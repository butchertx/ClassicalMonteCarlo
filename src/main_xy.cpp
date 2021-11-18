/**
Run a classical Monte Carlo simulation for n-rotor models with arbitrary interactions and arbitrary dimensions
**/

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <fstream>
#include "MemTimeTester.h"
#include <thread>
#include <chrono>
#include <mpi.h>

//#include "MCParams.h"
#include "cmctype.h"
#include "Lattice.h"
#include "MCState.h"
#include "RandomEngine.h"
#include "MCRun.h"
#include "MCResults.h"
#include "Update.h"


void check_gpu();

void write_results_mpi(int id, int p, cmctype::MCParams params, std::vector<double> process_results, std::string filename);

int main(int argc, char* argv[]) {
    //
    //  define global variables
    //
    int p = 1;//number of processes
    int id = 0;//ID of this process
    int seed = -1;//RNG seed for each process
    std::stringstream outstring;
//    int device_count = 0;//number of GPUs
//    int proc_per_device;//number of threads per GPU
//	bool gpu = true;//gpu, yay or nay, determined by command line argument
    //bool mpi_use = true;//mpi, yay or nay, determined by command line argument
//    cudaDeviceProp gpu_stats;
    

    /*
        if(argc >= 4){
            std::string arg1(argv[3]);
            if(arg1.compare("nompi") == 0){
                mpi_use = false;
            }
            else if(arg1.compare("nogpu") == 0){
                gpu = false;
            }
            if(argc >= 5){
                std::string arg2(argv[4]);
                if(arg2.compare("nompi") == 0){
                    mpi_use = false;
                }
                else if(arg2.compare("nogpu") == 0){
                    gpu = false;
                }
            }
        }
    */

    //std::cout << argc << "\n";
    //for (int i = 0; i < argc; ++i) {
    //    std::cout << "Arg " << i << ": " << argv[i] << "\n";
    //}
    bool mpi_use = false;
    for (int a = 0; a < argc; ++a) {
        std::string arg = argv[a];
        if (arg == "-mpi") {
            mpi_use = true;
        }
    }
    if (mpi_use) {
        //
        //  Initialize MPI.
        //
        MPI_Init(&argc, &argv);
        //
        //  Get the number of processes.
        //
        MPI_Comm_size(MPI_COMM_WORLD, &p);
        //
        //  Get the ID of this process.
        //
        MPI_Comm_rank(MPI_COMM_WORLD, &id);

        if (id == 0) {
            std::cout << "Using MPI with " << p << " processes\n";
        }
        std::cout << "Process " << id << " online\n";
    }
    else {
        p = 1;
        id = 0;
        std::cout << "Not using MPI\n";
    }

    //	if(id == 0){
    //        std::cout << "Using GPU, yes or no: " << (gpu ? "yes" : "no") << "\n";
    //    }


    //  Check GPU properties
    check_gpu();

    //
    //  Read and set up params for each process
    //
    cmctype::MCParams params;
    std::string infile(argv[argc-1]);
    params = cmctype::read_params(infile);
    params.set_parallel_params(id);
    params.print();
    std::vector<double> Jz1, beta;
    std::vector<double> mag2_results, mag_results, mag_abs_results;
    Jz1 = params.parallel.parallel_entries[1].values;

    //
    //  Start random number generators and set up lattice/state variables
    //  
    MemTimeTester timer;
    timer.flag_start_time("full simulation");
    Lattice lattice(Lattice_type_from_string(params.lattice.latticetype), params.lattice.L, vec3<int>(1, 1, 1));
    RandomEngine random = RandomEngine(params.markov.seed, lattice.get_N(), 1);
    assert(params.lattice.spintype.compare("XY") == 0);
    XYLattice state(lattice);
    Model_ANNNXY model(lattice, params.model);
    XYUpdate updater(&random, &model);
    MCRun<spin2> runner(&random, params, &state, &model, &updater);

    for (int j = 0; j < Jz1.size(); ++j) {
        if (params.parallel.parallel_entries[1].name == "Jz1") {
            params.model.interactions[1].strength = Jz1[j];
        }
        else {
            params.model.interactions[3].strength = Jz1[j];
        }
        params.model.T = 1 / params.model.beta;
        runner.reset_params(lattice, params);
        runner.reset_results();
        state.randomize(&random);
        runner.run();

        std::vector<double> mag2 = runner.get_results().get("m2");
        std::vector<double> mag = runner.get_results().get("m");
        std::vector<double> mag_abs = runner.get_results().get("|m|");
        std::vector<double> corr = runner.get_results().get_function_average("corr");
        mag2_results.push_back(0.0);
        mag_results.push_back(0.0);
        mag_abs_results.push_back(0.0);
        for (int i = 0; i < mag2.size(); ++i) {
            mag2_results[j] += mag2[i] / mag2.size();
            mag_results[j] += mag[i] / mag.size();
            mag_abs_results[j] += mag_abs[i] / mag_abs.size();
        }

        std::ofstream file;
        std::stringstream filename;
        if (params.parallel.parallel_entries[1].name == "Jz1") {
            filename << "corr_j" << j << "_beta" << id << ".csv";
        }
        else {
            filename << "corr_h" << j << "_beta" << id << ".csv";
        }
        file.open(filename.str());
        file << vec2str(corr);
    }

    if (mpi_use) {
        write_results_mpi(id, p, params, mag2_results, "mag2.csv");
        write_results_mpi(id, p, params, mag_results, "mag.csv");
        write_results_mpi(id, p, params, mag_abs_results, "mag_abs.csv");
    }


    timer.flag_end_time("full simulation");
    timer.print_timers();
    

    ////
    ////  End the MPI process
    ////
    if (mpi_use) {
        MPI_Finalize();
    }

    return 0;
}

void write_results_mpi(int id, int p, cmctype::MCParams params, std::vector<double> process_results, std::string filename) {
    if (id == 0) {
        std::vector<std::vector<double>> full_results;
        full_results.push_back(process_results);
        MPI_Status status;
        for (int i = 1; i < p; ++i) {
            //std::cout << "Receiving...#" << i << "\n";
            MPI_Recv(&process_results[0], process_results.size(), MPI_DOUBLE, i, 777, MPI_COMM_WORLD, &status);
            full_results.push_back(process_results);
        }
        std::vector<double> beta = params.parallel.parallel_entries[0].values;
        std::vector<double> other = params.parallel.parallel_entries[1].values;

        std::ofstream file;
        file.open(filename);

        if (params.parallel.parallel_entries[1].name == "Jz1") {
            file << "beta\\Jz1," << vec2str(other) << "\n";
        }
        else {
            file << "beta\\h," << vec2str(other) << "\n";
        }

        for (int b = 0; b < full_results.size(); ++b) {
            file << beta[b] << ",";
            for (int j = 0; j < other.size() - 1; ++j) {
                file << full_results[b][j] << ",";
            }
            file << full_results[b][other.size() - 1] << "\n";
        }
    }
    else {
        //std::cout << "Sending...\n";
        MPI_Send(&process_results[0], process_results.size(), MPI_DOUBLE, 0, 777, MPI_COMM_WORLD);
    }
}

void check_gpu() {
    /*
    if(gpu){
        cudaGetDeviceCount(&device_count);
        proc_per_device = p / device_count;
        if(id/proc_per_device < device_count){
            cudaSetDevice(id/proc_per_device);
        }
        else{
            cudaSetDevice(0);
        }
        if(id == 0){
            for(int i = 0; i < p; ++i){
                outstring << "Thread " << i << " set to Device " << i/proc_per_device << "\n";
            }
        }
        for(int d = 0; d < device_count; ++d){
            if(id == d * proc_per_device){
                cudaGetDeviceProperties(&gpu_stats, d);
                outstring << "GPU properties:\nDevice: " << d << "\nName: " << gpu_stats.name << "\nThreads per block: " << gpu_stats.maxThreadsPerBlock
                    << "\nThread dimensions: {" << gpu_stats.maxThreadsDim[0] << "," << gpu_stats.maxThreadsDim[1] << "," << gpu_stats.maxThreadsDim[2]
                    << "}\nMax grid size: {" << gpu_stats.maxGridSize[0] << "," << gpu_stats.maxGridSize[1] << "," << gpu_stats.maxGridSize[2] << "}\n\n";
                std::cout << outstring.str();
                outstring.str("");
                show_memory();
            }
        }
    }
*/
}