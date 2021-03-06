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

//#include "MCParams.h"
#include "cmctype.h"
#include "Lattice.h"
#include "MCState.h"
#include "RandomEngine.h"
#include "MCRun.h"
#include "MCResults.h"
#include "Update.h"


int main(int argc, char* argv[]) {

    int p = 1;//number of processes
    int id = 0;//ID of this process
    int seed = -1;//RNG seed for each process
    std::stringstream outstring;

    //
    //  Read and set up params for each process
    //
    cmctype::MCParams params;
    std::string infile(argv[argc-1]);
    params = cmctype::read_params(infile);
    std::vector<double> Jz1, beta;
    std::vector<std::vector<double>> mag2_results;
    beta = params.parallel.parallel_entries[0].values;
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
        mag2_results.push_back({});
        for (int b = 0; b < beta.size(); ++b) {
            params.model.interactions[1].strength = Jz1[j];
            params.model.beta = beta[b];
            params.model.T = 1 / beta[b];
            params.print();
            runner.reset_params(lattice, params);
            runner.reset_results();
            state.randomize(&random);
            runner.run();

            std::vector<double> mag2 = runner.get_results().get("mq");
            std::vector<double> corr = runner.get_results().get_function_average("corr");
            mag2_results[j].push_back(0.0);
            for (int i = 0; i < mag2.size(); ++i) {
                mag2_results[j][b] += mag2[i] / mag2.size();
            }

            std::ofstream file;
            std::stringstream filename;
            filename << "corr_j" << j << "_beta" << b << ".csv";
            file.open(filename.str());
            file << vec2str(corr);
        }
    }

    std::ofstream file;
    std::stringstream filename;
    filename << "mag2.csv";
    file.open(filename.str());
    file << "beta\\Jz1," << vec2str(Jz1) << "\n";
    for (int b = 0; b < mag2_results[0].size(); ++b) {
        file << beta[b] << ",";
        for (int j = 0; j < Jz1.size() - 1; ++j) {
            file << mag2_results[j][b] << ",";
        }
        file << mag2_results[Jz1.size() - 1][b] << "\n";
    }

    timer.flag_end_time("full simulation");
    timer.print_timers();

    return 0;
}