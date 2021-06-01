#pragma once

#include "json.hpp"
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

using json = nlohmann::json;

struct LatticeParams {

	std::string type;
	std::vector<int> L;
	int D;
	std::string pbc_type;

};

void to_json(json& j, const LatticeParams& p);

void from_json(const json& j, LatticeParams& p);

struct Interaction {
	
	std::string name;
	int distance;
	//std::vector<double> strength;
	double strength;

};

void to_json(json& j, const Interaction& p);

void from_json(const json& j, Interaction& p);

struct ModelParams {

	std::string interaction_type;
	std::vector<Interaction> interactions;
	double beta;
	double T;

};

void to_json(json& j, const ModelParams& p);

void from_json(const json& j, ModelParams& p);

struct ParallelEntry {

	std::string name;
	std::vector<double> values;
	bool mpi;

};

void to_json(json& j, const ParallelEntry& p);

void from_json(const json& j, ParallelEntry& p);

struct ParallelParams {

	std::vector<ParallelEntry> parallel_entries;

};

void to_json(json& j, const ParallelParams& p);

void from_json(const json& j, ParallelParams& p);

struct MarkovParams {

	int seed;
	int metropolis_steps;
	int overrelaxation_steps;
	int cluster_steps;
	int measures;
	int throwaway;
	int measures_per_ptemp;

};

void to_json(json& j, const MarkovParams& p);

void from_json(const json& j, MarkovParams& p);

class MCParams {

public:

	LatticeParams lattice;
	ModelParams model;
	ParallelParams parallel;
	MarkovParams markov;
	std::vector<std::string> observables;
	std::vector<std::string> observable_functions;
	std::vector<std::string> outputs;

	void print() {
		std::cout << "Lattice Params\n";
		std::cout << "Type: " << lattice.type << "\n";
		std::cout << "L: [";
		for (int i = 0; i < lattice.D-1; ++i) {
			std::cout << lattice.L[i] << ", ";
		}
		std::cout << lattice.L[lattice.L.size() - 1] << "]\n";
		std::cout << "Dimension: " << lattice.D << "\n";
		std::cout << "PBC type: " << lattice.pbc_type << "\n\n";

		std::cout << "Model Params\n";
		std::cout << "Interaction type: " << model.interaction_type << "\n";
		if (model.interaction_type.compare("local") == 0) {
			std::cout << "Interactions:\n";
		}
	}

	void set_parallel_params(int id);

};

MCParams read_params(json j);

MCParams read_params(std::string infile_name);