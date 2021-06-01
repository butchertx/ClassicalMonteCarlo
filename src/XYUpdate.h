#pragma once
#define _USE_MATH_DEFINES
#include <vector>
#include "RandomEngine.h"
#include "MCState.h"
#include <cmath>
#include <exception>



class XYUpdate {

protected:

	RandomEngine* random;
	XYUpdate(RandomEngine* rand_) : random(rand_) {};

public:

	virtual bool metropolis_step(XYLattice3D* lattice) = 0;

	virtual void overrelaxation_step(XYLattice3D* lattice) = 0;

	virtual void local_cluster_step(XYLattice3D* lattice) = 0;

};

class XYUpdate_3D_J0Jz1Jz2 : public XYUpdate {

	double K0;
	double Kz1;
	double Kz2;
	std::vector<std::vector<int>> neighbor_list = { {1, 0, 0}, 
													{-1, 0, 0}, 
													{0, 1, 0}, 
													{0, -1, 0}, 
													{0, 0, 1}, 
													{0, 0, -1}, 
													{0, 0, 2}, 
													{0, 0, -2} };
	std::vector<double> neighbor_couplings;

public:

	XYUpdate_3D_J0Jz1Jz2(RandomEngine* rand_, double K0_, double Kz1_, double Kz2_)
		: XYUpdate(rand_), K0(K0_), Kz1(Kz1_), Kz2(Kz2_) {
		neighbor_couplings = { K0, K0, K0, K0, Kz1, Kz1, Kz2, Kz2 };
		std::cout << "XYUpdate 3D J0Jz1Jz2 with K0 = " << K0 << ", Kz1 = " << Kz1 << ", Kz2 = " << Kz2 << "\n";
	}

	bool metropolis_step(XYLattice3D* lattice) {
		//select random spin and mirror axis
		int site0 = random->get_rand_site();
		spin2 seed = lattice->get(site0);
		double angle = random->get_rand_prob() * TAU;
		spin2 R = spin2(cos(angle), sin(angle));
		double seedR = seed.dot(R);

		//calculate energy change
		spin2 neigh;
		double action = 0.0;
		for (int n = 0; n < neighbor_list.size(); ++n) {
			neigh = lattice->get(site0, neighbor_list[n]);
			action += 2.0 * neighbor_couplings[n] * seedR * neigh.dot(R);
		}

		//accept or reject
		if (random->get_rand_prob() < exp(action)) {
			//should be able to refactor as a one-liner, maybe by a ternary operator
				lattice->set(site0, spin2(seed.sx - 2.0 * seedR * R.sx, seed.sy - 2.0 * seedR * R.sy));
			return true;
		}
		else {
			return false;
		}
	}

	void overrelaxation_step(XYLattice3D* lattice) {

		spin2 mean_field = spin2(0.0, 0.0);
		spin2 neigh, seed;
		double s_dot_mf;

		//iterate over sites
		for (int site0 = 0; site0 < lattice->size(); ++site0) {

			seed = lattice->get(site0);
			mean_field = spin2(0.0, 0.0);

			//Calculate mean field around site0
			for (int n = 0; n < neighbor_list.size(); ++n) {
				neigh = lattice->get(site0, neighbor_list[n]);
				mean_field += (neigh * neighbor_couplings[n]);
			}

			//rotate without changing energy
			mean_field.normalize();
			s_dot_mf = mean_field.dot(seed);
			lattice->set(site0, spin2(seed.sx - 2.0 * s_dot_mf * mean_field.sx, seed.sy - 2.0 * s_dot_mf * mean_field.sy));
		}
	}

	void local_cluster_step(XYLattice3D* lattice) {
		//create buffer and start cluster
		std::vector<int> buffer, flip_list;
		std::vector<bool> cluster(lattice->size(), false);
		buffer.push_back(random->get_rand_site());
		double angle = random->get_rand_prob() * TAU;
		spin2 R = spin2(cos(angle), sin(angle));
		
		int site0 = buffer.back();
		cluster[site0] = true;
		flip_list.push_back(site0);
		spin2 seed, neigh;
		double seedR;
		while (buffer.size() > 0) {
			site0 = buffer.back();
			buffer.pop_back();

			seed = lattice->get(site0);
			seedR = seed.dot(R);
			for (int n = 0; n < neighbor_list.size(); ++n) {
				neigh = lattice->get(site0, neighbor_list[n]);
				if (!cluster[neigh.site]
					&& random->get_rand_prob() < 1 - exp(2.0 * neighbor_couplings[n] * seedR * neigh.dot(R))
					){//&& neighbor_couplings[n] * seed.dot(neigh) < 0) {
					buffer.push_back(neigh.site);
					flip_list.push_back(neigh.site);
					cluster[neigh.site] = true;
				}
			}
		}

		for (int i = 0; i < flip_list.size(); ++i) {
			seed = lattice->get(flip_list[i]);
			seedR = seed.dot(R);
			lattice->set(flip_list[i], spin2(seed.sx - 2.0 * seedR * R.sx, seed.sy - 2.0 * seedR * R.sy));
		}


	};


};