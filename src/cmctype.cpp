#include "cmctype.h"

cmctype::MCParams cmctype::read_params(json j) {
	cmctype::MCParams params;

	params.lattice = j["lattice"].get<cmctype::LatticeParams>();
	params.model = j["model"].get<cmctype::ModelParams>();
	params.parallel = j["parallel"].get<cmctype::ParallelParams>();
	params.markov = j["markov"].get<cmctype::MarkovParams>();
	params.observables = j["observables"].get<std::vector<std::string>>();
	params.observable_functions = j["observable_functions"].get<std::vector<std::string>>();
	params.outputs = j["outputs"].get<std::vector<std::string>>();

	return params;
}

cmctype::MCParams cmctype::read_params(std::string infile_name) {
	std::ifstream i(infile_name);
	json j;
	i >> j;
	return cmctype::read_params(j);
}

void cmctype::MCParams::print(){
	std::cout << "Lattice Params\n";
	std::cout << "Spin Type: " << lattice.spintype << "\n";
	std::cout << "Lattice Type: " << lattice.latticetype << "\n";
	std::cout << "L: [";
	for (int i = 0; i < lattice.L.size() - 1; ++i) {
		std::cout << lattice.L[i] << ", ";
	}
	std::cout << lattice.L[lattice.L.size() - 1] << "]\n";
	std::cout << "PBC type: " << lattice.pbc_type << "\n\n";

	std::cout << "Model Params\n";
	std::cout << "Model: " << model.model_name << "\n";
	std::cout << "Interactions:\n";
}

void cmctype::to_json(json& j, const cmctype::LatticeParams& p) {
	j = json{
		{"spintype", p.spintype},
		{"latticetype", p.latticetype},
		{"L", p.L},
		{"pbc_type", p.pbc_type}
	};
}

void cmctype::from_json(const json& j, cmctype::LatticeParams& p) {
	j.at("spintype").get_to(p.spintype);
	j.at("latticetype").get_to(p.latticetype);
	j.at("L").get_to<std::vector<int>>(p.L);
	j.at("pbc_type").get_to(p.pbc_type);
}

void cmctype::to_json(json& j, const cmctype::Interaction& p) {
	j = json{
		{"name", p.name},
		{"strength", p.strength}
	};
}

void cmctype::from_json(const json& j, cmctype::Interaction& p) {
	j.at("name").get_to(p.name);
	j.at("strength").get_to(p.strength);
}

void cmctype::to_json(json& j, const cmctype::ModelParams& p) {
	j = json{
		{"model_name", p.model_name},
		{"interactions", p.interactions},
		{"beta", p.beta}
	};
}

void cmctype::from_json(const json& j, cmctype::ModelParams& p) {
	j.at("model_name").get_to(p.model_name);
	j.at("interactions").get_to<std::vector<cmctype::Interaction>>(p.interactions);
	j.at("beta").get_to(p.beta);
	p.T = 1.0 / p.beta;
}

void cmctype::to_json(json& j, const cmctype::ParallelEntry& p) {
	j = json{
		{"name", p.name},
		{"values", p.values},
		{"MPI", p.mpi}
	};
}

void cmctype::from_json(const json& j, cmctype::ParallelEntry& p) {
	j.at("name").get_to(p.name);
	j.at("values").get_to<std::vector<double>>(p.values);
	j.at("MPI").get_to(p.mpi);
}

void cmctype::to_json(json& j, const cmctype::ParallelParams& p) {
	j = json{
		{"parallel_entries", p.parallel_entries}
	};
}

void cmctype::from_json(const json& j, cmctype::ParallelParams& p) {
	j.at("parallel_entries").get_to<std::vector<cmctype::ParallelEntry>>(p.parallel_entries);
}

void cmctype::to_json(json& j, const cmctype::MarkovParams& p) {
	j = json{
		{"seed", p.seed},
		{"metropolis_steps", p.metropolis_steps},
		{"overrelaxation_steps", p.overrelaxation_steps},
		{"cluster_steps", p.cluster_steps},
		{"measures", p.measures},
		{"throwaway", p.throwaway},
		{"measures_per_ptemp", p.measures_per_ptemp}
	};
}

void cmctype::from_json(const json& j, cmctype::MarkovParams& p) {
	j.at("seed").get_to(p.seed);
	j.at("metropolis_steps").get_to(p.metropolis_steps);
	j.at("overrelaxation_steps").get_to(p.overrelaxation_steps);
	j.at("cluster_steps").get_to(p.cluster_steps);
	j.at("measures").get_to(p.measures);
	j.at("throwaway").get_to(p.throwaway);
	j.at("measures_per_ptemp").get_to(p.measures_per_ptemp);
}

void cmctype::MCParams::set_parallel_params(int id) {
	std::cout << "Parallel params not implemented yet\n";
}
