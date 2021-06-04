#pragma once
#include <unordered_map>
#include <string>
#include <assert.h>
#include "MCResults.h"
#include "RandomEngine.h"
#include "MCParams.h"
#include "MCState.h"
#include "XYUpdate.h"
#include <typeinfo>
#include <vector>
#ifndef _USE_MATH_DEFINES
	#define _USE_MATH_DEFINES
#endif


class MCRun {
	
	std::unordered_map<std::string, int> step_counts;
	int measure_target;
	int measures_per_ptemp;
	RandomEngine* rand_p;
	MCResults results;
	MCParams params;
	XYLattice3D* lattice;

	std::vector<double> Q = { 0.0, 0.0, 0.5 * TAU };

public:

	MCRun(RandomEngine* rand_, MCParams params_, XYLattice3D* lattice_): rand_p(rand_), params(params_), lattice(lattice_) {
		if (params.model.interactions[1].strength > 0) {
			Q = { 0.0, 0.0, 0.5 * TAU };
		}
		else {
			Q = { 0.0, 0.0, 0.0 };
		}
	}

	void run() {

		XYUpdate_3D_J0Jz1Jz2 updater(rand_p,
			params.model.beta * params.model.interactions[0].strength,
			params.model.beta * params.model.interactions[1].strength,
			params.model.beta * params.model.interactions[2].strength);
		for (int m = 0; m < params.markov.measures + params.markov.throwaway; ++m) {
			for (int i = 0; i < params.markov.overrelaxation_steps; ++i) {
				updater.overrelaxation_step(lattice);
			}
			for (int i = 0; i < params.markov.metropolis_steps; ++i) {
				updater.metropolis_step(lattice);
			}
			for (int i = 0; i < params.markov.cluster_steps; ++i) {
				updater.local_cluster_step(lattice);
			}
			if (m > params.markov.throwaway) {
				measure();
			}
		}

		for (int i = 0; i < lattice->size(); ++i) {
			assert(lattice->get(i).site == i);
		}
	}

	MCResults get_results() {
		return results;
	}

	void reset_params(MCParams params_) {
		params = params_;
		if (params.model.interactions[1].strength > 0) {
			Q = { 0.0, 0.0, 0.5 * TAU };
		}
		else {
			Q = { 0.0, 0.0 };
		}
	}

	void reset_results() {
		results = MCResults();
	}

	void measure() {

		for (auto obs : params.observables) {
			results.record(obs, calc_observable(obs));
		}
		for (auto obs : params.observable_functions) {
			results.record(obs, calc_observable_function(obs));
		}
	}

	double calc_observable(std::string obs_name) {
		if (obs_name.compare("m2") == 0) {
			return lattice->calc_mag2();
		}
		else if (obs_name.compare("mq") == 0) {
			return lattice->calc_chi_q(Q);
		}
		else {
			std::cout << obs_name << " not supported as observable\n";
			assert(1 == 0);
		}
	}

	std::vector<double> calc_observable_function(std::string obs_name) {
		if (obs_name.compare("corr") == 0) {
			return lattice->calc_correlation();
		}
		else {
			std::cout << obs_name << " not supported as observable\n";
			assert(1 == 0);
		}
	}

};