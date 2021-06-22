#pragma once
#include <assert.h>

#include "cmctype.h"
#include "Lattice.h"

class Model {

public:

	//index 1: coupling, index 2: origin site, index 3: neighbor sites
	std::vector<std::vector<std::vector<int>>> neighbor_list;
	//index 1: coupling strength
	std::vector<double> neighbor_couplings;

	vec3<double> Q = { 0.0, 0.0, 0.0 };

	Model() {};

};

class Model_ANNNXY : public Model {

public:

	Model_ANNNXY(Lattice& lattice_, cmctype::ModelParams params_) {
		assert(lattice_.get_lattice_type() == CUBIC);
		assert(params_.model_name.compare("ANNNXY") == 0);

		std::vector<vec3<double>> displacements = { vec3<double>(-1.0, 0.0, 0.0),
													vec3<double>(1.0, 0.0, 0.0),
													vec3<double>(0.0, -1.0, 0.0),
													vec3<double>(0.0, 1.0, 0.0),
													vec3<double>(0.0, 0.0, -1.0),
													vec3<double>(0.0, 0.0, 1.0),
													vec3<double>(0.0, 0.0, -2.0),
													vec3<double>(0.0, 0.0, 2.0) };

		neighbor_list = { {}, {}, {} };
		neighbor_couplings = { params_.beta * params_.get_interaction("J0"), params_.beta * params_.get_interaction("Jz1"), params_.beta * params_.get_interaction("Jz2") };

		for (int i = 0; i < lattice_.get_N(); ++i) {
			
			neighbor_list[0].push_back({});
			neighbor_list[0][i].push_back(lattice_.get_neighbor_label(i, displacements[0]));
			neighbor_list[0][i].push_back(lattice_.get_neighbor_label(i, displacements[1]));
			neighbor_list[0][i].push_back(lattice_.get_neighbor_label(i, displacements[2]));
			neighbor_list[0][i].push_back(lattice_.get_neighbor_label(i, displacements[3]));

			neighbor_list[1].push_back({});
			neighbor_list[1][i].push_back(lattice_.get_neighbor_label(i, displacements[4]));
			neighbor_list[1][i].push_back(lattice_.get_neighbor_label(i, displacements[5]));

			neighbor_list[2].push_back({});
			neighbor_list[2][i].push_back(lattice_.get_neighbor_label(i, displacements[6]));
			neighbor_list[2][i].push_back(lattice_.get_neighbor_label(i, displacements[7]));
		}
	}

};