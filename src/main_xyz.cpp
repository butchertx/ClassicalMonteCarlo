/**
Run a classical Monte Carlo simulation for n-rotor models with arbitrary interactions and arbitrary dimensions
**/

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <thread>
#include <chrono>

#include "MemTimeTester.h"
#include "cmctype.h"
#include "Lattice.h"
#include "MCState.h"
#include "RandomEngine.h"
#include "MCRun.h"
#include "MCResults.h"
#include "Update.h"
#include "Model.h"


int main(int argc, char* argv[]) {
    //
    //  define global variables
    //
    int p = 1;//number of processes
    int id = 0;//ID of this process
    int seed = -1;//RNG seed for each process
    std::stringstream outstring;

    //
    //  Read and set up params for each process
    //
    cmctype::MCParams params;
    std::string infile(argv[1]);
    params = cmctype::read_params(infile);
    params.set_parallel_params(id);
    params.print();
    std::vector<double> beta;
    std::vector<std::vector<double>> mag2_results;
    beta = params.parallel.parallel_entries[0].values;

    //
    //  Start random number generators and set up lattice/state variables
    //

    MemTimeTester timer;
    timer.flag_start_time("full simulation");
    Lattice lattice(Lattice_type_from_string(params.lattice.latticetype), params.lattice.L, vec3<int>(1, 1, 1));
    RandomEngine random = RandomEngine(params.markov.seed, lattice.get_N(), 1);
    assert(params.lattice.spintype.compare("XYZ") == 0);
    XYZLattice state(lattice);
    Model_CGT model(lattice, params.model);
    
    XYZUpdate updater(&random, &model);
    MCRun<spin3> runner(&random, params, &state, &model, &updater);

    mag2_results.push_back({});
    for (int b = 0; b < beta.size(); ++b) {
        params.model.beta = beta[b];
        params.model.T = 1 / beta[b];
        //runner.reset_params(params);
        model.set_interactions(lattice, params.model);
        runner.reset_results();
        state.randomize(&random);
        runner.run();

        std::vector<double> mag2 = runner.get_results().get("m2");
        std::vector<double> corr = runner.get_results().get_function_average("corr");
        mag2_results[0].push_back(0.0);
        for (int i = 0; i < mag2.size(); ++i) {
            mag2_results[0][b] += mag2[i] / mag2.size();
        }

        std::ofstream file;
        std::stringstream filename;
        filename << "corr_beta" << b << ".csv";
        file.open(filename.str());
        file << vec2str(corr);
    }

    std::ofstream file;
    file.open("mag2.csv");
    file << "beta, m^2\n";
    for (int b = 0; b < beta.size(); ++b) {
        file << beta[b] << ",";
        file << mag2_results[0][b] << "\n";
    }


    timer.flag_end_time("full simulation");
    timer.print_timers();

    return 0;
}