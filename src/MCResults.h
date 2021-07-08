#pragma once
#include <string>
#include <vector>
#include <sstream>
#include <unordered_map>


template <class T>
std::string vec2str(std::vector<T> vec) {
	std::stringstream ss;
	for (int i = 0; i < vec.size() - 1; ++i) {
		ss << vec[i] << ", ";
	}
	ss << vec.back();
	return ss.str();
}

bool isDirExist(const std::string& path);

bool makePath(const std::string& path);

class MCResults {

	std::unordered_map<std::string, std::vector<double>> measurements;
	std::unordered_map<std::string, std::vector<std::vector<double>>> function_measurements;

	std::unordered_map<std::string, std::vector<std::vector<double>>> function_bins;
	std::unordered_map<std::string, std::vector<std::vector<double>>> function_sq_bins;
	std::unordered_map<std::string, std::vector<int>> function_counts;

	bool keep_functions = false;
	int function_bin_size = 1;
	int measure_bin_size = 100;

public:

	MCResults() {
		makePath("./dump");

		/*
		//make path for write states and write process params
		
		if (params.alg.compare("write_states") == 0) {
		    char state_dump[100];
		    sprintf(state_dump, "dump/statedump%d", id);
		    makePath(state_dump);
		}
		char param_file_name[100];
		ofstream param_file_out;
		sprintf(param_file_name, "dump/params%d.dat", id);
		param_file_out.open(param_file_name);
		param_file_out << params.to_string();
		param_file_out.close();
		*/

	}

	void record(std::string obs, double measurement);

	void record(std::string obs, std::vector<double> measurement);

	void bin_function(std::string obs);

	std::vector<double> get(std::string obs);

	std::vector<double> get_function_average(std::string obs);

	std::vector<double> get_function_var(std::string obs);

	std::vector<double> autocorrelation(std::vector<double> measurements);

};