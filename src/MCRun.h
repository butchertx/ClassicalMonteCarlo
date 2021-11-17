#pragma once
#include <unordered_map>
#include <string>
#include <assert.h>
#include "cmctype.h"
#include "MCResults.h"
#include "RandomEngine.h"
#include "MCState.h"
#include "Update.h"
#include <typeinfo>
#include <vector>
#ifndef _USE_MATH_DEFINES
	#define _USE_MATH_DEFINES
#endif

template<class spintype>
class MCRun {
	
	std::unordered_map<std::string, int> step_counts;
	int measure_target;
	int measures_per_ptemp;
	RandomEngine* rand_p;
	MCResults results;
	cmctype::MCParams params;
	NRotorLattice<spintype>* state;
	Model* model;
	MCUpdate<spintype>* updater;

public:

	MCRun(RandomEngine* rand_, cmctype::MCParams params_, NRotorLattice<spintype>* state_, Model* model_, MCUpdate<spintype>* updater_) : rand_p(rand_), params(params_), state(state_), model(model_), updater(updater_){};

	void run() {

		for (int m = 0; m < params.markov.measures + params.markov.throwaway; ++m) {
			for (int i = 0; i < params.markov.overrelaxation_steps; ++i) {
				updater->overrelaxation_step(state);
			}
			for (int i = 0; i < params.markov.metropolis_steps; ++i) {
				updater->metropolis_step(state);
			}
			for (int i = 0; i < params.markov.cluster_steps; ++i) {
				updater->local_cluster_step(state);
			}
			if (m > params.markov.throwaway) {
				measure();
			}
		}

		for (int i = 0; i < state->size(); ++i) {
			assert(state->get(i).site == i);
		}
	}

	MCResults get_results() {
		return results;
	}

	void reset_params(Lattice& lattice_, cmctype::MCParams params_) {
		params = params_;
		model->reset_interactions(lattice_, params_.model);
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
			return state->calc_mag2();
		}
		else if (obs_name.compare("mq") == 0) {
			return state->calc_chi_q(model->Q);
		}
		else {
			std::cout << obs_name << " not supported as observable\n";
			assert(1 == 0);
		}
	}

	std::vector<double> calc_observable_function(std::string obs_name) {
		if (obs_name.compare("corr") == 0) {
			return state->calc_correlation();
		}
		else {
			std::cout << obs_name << " not supported as observable\n";
			assert(1 == 0);
		}
	}

};