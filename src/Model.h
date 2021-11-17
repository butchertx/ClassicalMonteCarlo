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

	Model(Lattice& lattice_, cmctype::ModelParams params_) {
		set_interactions(lattice_, params_);
	}

	void virtual set_interactions(Lattice& lattice_, cmctype::ModelParams params_) {};

	void reset_interactions(Lattice& lattice_, cmctype::ModelParams params_) {
		set_interactions(lattice_, params_);
	}

};

class Model_ANNNXY : public Model {

public:

	Model_ANNNXY(Lattice& lattice_, cmctype::ModelParams params_) : Model(lattice_, params_) {
		set_interactions(lattice_, params_);
	};

	void set_interactions(Lattice& lattice_, cmctype::ModelParams params_) {
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
	};

};

class Model_Heisenberg : public Model {

public:

	Model_Heisenberg(Lattice& lattice_, cmctype::ModelParams params_) : Model(lattice_, params_) {
		set_interactions(lattice_, params_);
	};

	void set_interactions(Lattice & lattice_, cmctype::ModelParams params_) {

		assert(lattice_.get_lattice_type() == CUBIC);
		assert(params_.model_name.compare("Heisenberg") == 0);

		std::vector<vec3<double>> displacements = { vec3<double>(-1.0, 0.0, 0.0),
													vec3<double>(1.0, 0.0, 0.0),
													vec3<double>(0.0, -1.0, 0.0),
													vec3<double>(0.0, 1.0, 0.0),
													vec3<double>(0.0, 0.0, -1.0),
													vec3<double>(0.0, 0.0, 1.0) };

		neighbor_list = { {} };
		neighbor_couplings = { params_.beta * params_.get_interaction("J1")};

		for (int i = 0; i < lattice_.get_N(); ++i) {

			neighbor_list[0].push_back({});
			neighbor_list[0][i].push_back(lattice_.get_neighbor_label(i, displacements[0]));
			neighbor_list[0][i].push_back(lattice_.get_neighbor_label(i, displacements[1]));
			neighbor_list[0][i].push_back(lattice_.get_neighbor_label(i, displacements[2]));
			neighbor_list[0][i].push_back(lattice_.get_neighbor_label(i, displacements[3]));
			neighbor_list[0][i].push_back(lattice_.get_neighbor_label(i, displacements[4]));
			neighbor_list[0][i].push_back(lattice_.get_neighbor_label(i, displacements[5]));
		}
	}

};

class Model_CGT : public Model {

public:

	Model_CGT(Lattice& lattice_, cmctype::ModelParams params_) : Model(lattice_, params_) {
		set_interactions(lattice_, params_);
	};

	void set_interactions(Lattice& lattice_, cmctype::ModelParams params_) {

		assert(lattice_.get_lattice_type() == CGT);
		assert(params_.model_name.compare("CGT") == 0);

		std::vector<vec3<double>> basis = lattice_.get_basis();
		vec3<vec3<double>> a = lattice_.get_lattice_vectors();

		std::vector<vec3<double>> displacements = { basis[1], //J1 (A origin)
													basis[1] - a.y, //J1 (A origin)
													basis[1] + a.x - a.y, //J1 (A origin)
													basis[1] * (-1.0), //J1 (B origin)
													a.y - basis[1], //J1 (B origin)
													a.y - basis[1] - a.x, //J1 (B origin)

													a.x, //J2
													a.y, //J2
													a.y - a.x, //J2
													a.x * (-1.0), //J2
													a.y * (-1.0), //J2
													a.x - a.y, //J2

													a.x + basis[1], //J3 (A origin)
													basis[1] - a.x, //J3 (A origin)
													a.x - (a.y * (2.0)) + basis[1], //J3 (A origin)
													a.x * (-1.0) - basis[1], //J3 (B origin)
													a.x - basis[1], //J3 (B origin)
													(a.y * (2.0)) - a.x - basis[1], //J3 (B origin)

													basis[1] - a.z, //Jz1 (A origin)
													a.z - basis[1], //Jz1 (B origin)

													a.z + a.x - a.y + basis[1], //Jz2 (A origin)
													a.z - a.y + basis[1], //Jz2 (A origin)
													a.z - a.y * (2.0) + a.x + basis[1], //Jz2 (A origin)
													a.y - a.z - a.x - basis[1], //Jz2 (B origin)
													a.y - a.z - basis[1], //Jz2 (B origin)
													a.y * (2.0) - a.z - a.x - basis[1], //Jz2 (B origin)

													a.z, //Jz3
													a.z - a.y, //Jz3
													a.z - a.y + a.x, //Jz3
													a.z * (-1.0), //Jz3
													a.y - a.z, //Jz3
													a.y - a.z - a.x, //Jz3
		};

		std::vector<int> displacement_coupling_indices = { 0, 0, 0, 0, 0, 0, //J1
														1, 1, 1, 1, 1, 1, //J2
														2, 2, 2, 2, 2, 2, //J3
														3, 3, //Jz1
														4, 4, 4, 4, 4, 4, //Jz2
														5, 5, 5, 5, 5, 5 }; //Jz3

		assert(displacement_coupling_indices.size() == displacements.size());
		neighbor_list = { {}, {}, {}, {}, {}, {} };
		neighbor_couplings = { params_.beta * params_.get_interaction("J1"), params_.beta * params_.get_interaction("J2"), params_.beta * params_.get_interaction("J3"),
								params_.beta * params_.get_interaction("Jz1"), params_.beta * params_.get_interaction("Jz2"), params_.beta * params_.get_interaction("Jz3") };


		
		for (int i = 0; i < lattice_.get_N(); ++i) {
			for (int n = 0; n < neighbor_couplings.size(); ++n) {
				neighbor_list[n].push_back({});
			}
			for (int n = 0; n < displacement_coupling_indices.size(); ++n) {
				if (lattice_.check_shift(i, displacements[n])) {
					neighbor_list[displacement_coupling_indices[n]][i].push_back(lattice_.get_neighbor_label(i, displacements[n]));
				}
			}
		}
	}

};