#pragma once
#define _USE_MATH_DEFINES
#include <vector>
#include "RandomEngine.h"
#include "MCState.h"
#include <cmath>
#include <exception>
#include "Model.h"

template <class spintype>
class MCUpdate {

protected:

	RandomEngine* random;
	Model* model;
	MCUpdate(RandomEngine* rand_, Model* model_) : random(rand_), model(model_) {};

public:

	virtual bool metropolis_step(NRotorLattice<spintype>* state) = 0;

	virtual void overrelaxation_step(NRotorLattice<spintype>* state) = 0;

	virtual void local_cluster_step(NRotorLattice<spintype>* state) = 0;

};

class XYUpdate : public MCUpdate<spin2> {

public:

	XYUpdate(RandomEngine* rand_, Model* model_) : MCUpdate<spin2>(rand_, model_) {};

	bool metropolis_step(NRotorLattice<spin2>* state) {
		//select random spin and mirror axis
		int site0 = random->get_rand_site();
		spin2 seed = state->get(site0);
		double angle = random->get_rand_prob() * TAU;
		cmctype::vec2<double> R(cos(angle), sin(angle));
		double seedR = seed.dot(R);

		//calculate energy change
		spin2 neigh;
		double action = 0.0;
		for (int n = 0; n < model->neighbor_couplings.size(); ++n) {
			for (int other : model->neighbor_list[n][site0]) {
				neigh = state->get(other);
				action += 2.0 * model->neighbor_couplings[n] * seedR * neigh.dot(R);
			}
		}

		//accept or reject
		if (random->get_rand_prob() < exp(action)) {
			state->set(site0, seed - (R*(2.0 * seedR)));
			return true;
		}
		else {
			return false;
		}
	}

	void overrelaxation_step(NRotorLattice<spin2>* state) {

		spin2 mean_field(0.0, 0.0);
		spin2 neigh, seed;
		double s_dot_mf;
		int site0 = random->get_rand_site();
		seed = state->get(site0);

		//Calculate mean field around site0
		for (int n = 0; n < model->neighbor_couplings.size(); ++n) {
			for (int other : model->neighbor_list[n][site0]) {
				neigh = state->get(other);
				mean_field += (neigh * model->neighbor_couplings[n]);
			}
		}

		//rotate without changing energy
		mean_field.normalize();
		s_dot_mf = mean_field * seed;
		state->set(site0, seed - (mean_field * (2.0 * s_dot_mf)));

	}

	void local_cluster_step(NRotorLattice<spin2>* state) {
		//create buffer and start cluster
		std::vector<int> buffer, flip_list;
		std::vector<bool> cluster(state->size(), false);
		buffer.push_back(random->get_rand_site());
		double angle = random->get_rand_prob() * TAU;
		cmctype::vec2<double> R(cos(angle), sin(angle));
		
		int site0 = buffer.back();
		cluster[site0] = true;
		flip_list.push_back(site0);
		spin2 seed, neigh;
		double seedR;
		while (buffer.size() > 0) {
			site0 = buffer.back();
			buffer.pop_back();

			seed = state->get(site0);
			seedR = seed * R;
			for (int n = 0; n < model->neighbor_couplings.size(); ++n) {
				for (int other : model->neighbor_list[n][site0]) {
					neigh = state->get(other);
					if (!cluster[neigh.site]
						&& random->get_rand_prob() < 1 - exp(2.0 * model->neighbor_couplings[n] * seedR * neigh.dot(R))
						&& (seed * neigh * (model->neighbor_couplings[n]))  < 0) {
						buffer.push_back(neigh.site);
						flip_list.push_back(neigh.site);
						cluster[neigh.site] = true;
					}
				}
			}
		}

		for (int i = 0; i < flip_list.size(); ++i) {
			seed = state->get(flip_list[i]);
			seedR = seed * R;
			state->set(flip_list[i], seed - R * (2.0 * seedR));
		}


	};


};

class XYZUpdate : public MCUpdate<spin3> {

public:

	XYZUpdate(RandomEngine* rand_, Model* model_) : MCUpdate<spin3>(rand_, model_) {};

	bool metropolis_step(NRotorLattice<spin3>* state) {
		//select random spin and mirror axis
		int site0 = random->get_rand_site();
		spin3 seed = state->get(site0);
		double theta = random->get_rand_prob() * 0.5 * TAU;
		double phi = random->get_rand_prob() * TAU;
		cmctype::vec3<double> R(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
		double seedR = seed.dot(R);

		//calculate energy change
		spin3 neigh;
		double action = 0.0;
		for (int n = 0; n < model->neighbor_couplings.size(); ++n) {
			for (int other : model->neighbor_list[n][site0]) {
				neigh = state->get(other);
				action += 2.0 * model->neighbor_couplings[n] * seedR * neigh.dot(R);
			}
		}

		//accept or reject
		if (random->get_rand_prob() < exp(action)) {
			state->set(site0, seed - (R * (2.0 * seedR)));
			return true;
		}
		else {
			return false;
		}
	}

	void overrelaxation_step(NRotorLattice<spin3>* state) {

		spin3 mean_field(0.0, 0.0, 0.0);
		spin3 neigh, seed;
		double s_dot_mf;
		int site0 = random->get_rand_site();
		seed = state->get(site0);

		//Calculate mean field around site0
		for (int n = 0; n < model->neighbor_couplings.size(); ++n) {
			for (int other : model->neighbor_list[n][site0]) {
				neigh = state->get(other);
				mean_field += (neigh * model->neighbor_couplings[n]);
			}
		}

		//rotate without changing energy
		mean_field.normalize();
		s_dot_mf = mean_field * seed;
		state->set(site0, seed - (mean_field * (2.0 * s_dot_mf)));

	}

	void local_cluster_step(NRotorLattice<spin3>* state) {
		//create buffer and start cluster
		std::vector<int> buffer, flip_list;
		std::vector<bool> cluster(state->size(), false);
		buffer.push_back(random->get_rand_site());
		double theta = random->get_rand_prob() * 0.5 * TAU;
		double phi = random->get_rand_prob() * TAU;
		cmctype::vec3<double> R(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));

		int site0 = buffer.back();
		cluster[site0] = true;
		flip_list.push_back(site0);
		spin3 seed, neigh;
		double seedR;
		while (buffer.size() > 0) {
			site0 = buffer.back();
			buffer.pop_back();

			seed = state->get(site0);
			seedR = seed * R;
			for (int n = 0; n < model->neighbor_couplings.size(); ++n) {
				for (int other : model->neighbor_list[n][site0]) {
					neigh = state->get(other);
					if (!cluster[neigh.site]
						&& random->get_rand_prob() < 1 - exp(2.0 * model->neighbor_couplings[n] * seedR * neigh.dot(R))
						&& (seed * neigh * (model->neighbor_couplings[n])) < 0) {
						buffer.push_back(neigh.site);
						flip_list.push_back(neigh.site);
						cluster[neigh.site] = true;
					}
				}
			}
		}

		for (int i = 0; i < flip_list.size(); ++i) {
			seed = state->get(flip_list[i]);
			seedR = seed * R;
			state->set(flip_list[i], seed - R * (2.0 * seedR));
		}


	};


};