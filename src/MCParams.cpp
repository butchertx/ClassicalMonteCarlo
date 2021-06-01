#include "MCParams.h"
#include <string>

MCParams read_params(json j) {
	MCParams params;

	params.lattice = j["lattice"].get<LatticeParams>();
	params.model = j["model"].get<ModelParams>();
	params.parallel = j["parallel"].get<ParallelParams>();
	params.markov = j["markov"].get<MarkovParams>();
	params.observables = j["observables"].get<std::vector<std::string>>();
	params.observable_functions = j["observable_functions"].get<std::vector<std::string>>();
	params.outputs = j["outputs"].get<std::vector<std::string>>();

	return params;
}

MCParams read_params(std::string infile_name) {
	std::ifstream i(infile_name);
	json j;
	i >> j;
	return read_params(j);
}

void to_json(json& j, const LatticeParams& p) {
	j = json{
			{"type", p.type},
			{"L", p.L},
			{"dimension", p.D},
			{"pbc_type", p.pbc_type}
	};
}

void from_json(const json& j, LatticeParams& p) {
	j.at("type").get_to(p.type);
	j.at("L").get_to<std::vector<int>>(p.L);
	j.at("dimension").get_to(p.D);
	j.at("pbc_type").get_to(p.pbc_type);
}

void to_json(json& j, const Interaction& p) {
	j = json{
			{"name", p.name},
			{"distance", p.distance},
			{"strength", p.strength}
	};
}

void from_json(const json& j, Interaction& p) {
	j.at("name").get_to(p.name);
	j.at("distance").get_to(p.distance);
	//j.at("strength").get_to<std::vector<double>>(p.strength);
	j.at("strength").get_to(p.strength);
}

void to_json(json& j, const ModelParams& p) {
	j = json{
			{"interaction_type", p.interaction_type},
			{"interactions", p.interactions},
			{"beta", p.beta}
	};
}

void from_json(const json& j, ModelParams& p) {
	j.at("interaction_type").get_to(p.interaction_type);
	if (p.interaction_type.compare("local") == 0) {
		j.at("interactions").get_to<std::vector<Interaction>>(p.interactions);
	}
	j.at("beta").get_to(p.beta);
	p.T = 1.0/p.beta;
}

void to_json(json& j, const ParallelEntry& p) {
	j = json{
		{"name", p.name},
		{"values", p.values},
		{"MPI", p.mpi}
	};
}

void from_json(const json& j, ParallelEntry& p) {
	j.at("name").get_to(p.name);
	j.at("values").get_to<std::vector<double>>(p.values);
	j.at("MPI").get_to(p.mpi);
}

void to_json(json& j, const ParallelParams& p) {
	j = json{
		{"parallel_entries", p.parallel_entries}
	};
}

void from_json(const json& j, ParallelParams& p) {
	j.at("parallel_entries").get_to<std::vector<ParallelEntry>>(p.parallel_entries);
}

void to_json(json& j, const MarkovParams& p) {
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

void from_json(const json& j, MarkovParams& p) {
	j.at("seed").get_to(p.seed);
	j.at("metropolis_steps").get_to(p.metropolis_steps);
	j.at("overrelaxation_steps").get_to(p.overrelaxation_steps);
	j.at("cluster_steps").get_to(p.cluster_steps);
	j.at("measures").get_to(p.measures);
	j.at("throwaway").get_to(p.throwaway);
	j.at("measures_per_ptemp").get_to(p.measures_per_ptemp);
}

void MCParams::set_parallel_params(int id) {
	std::cout << "Parallel params not implemented yet\n";
}
