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

	void record(std::string obs, double measurement) {
		measurements[obs].push_back(measurement);
	}

	void record(std::string obs, std::vector<double> measurement) {
		function_measurements[obs].push_back(measurement);
	}

	std::vector<double> get(std::string obs) {
		return measurements[obs];
	}

	std::vector<double> get_function_average(std::string obs) {
		std::vector<std::vector<double>> results = function_measurements[obs];
		std::vector<double> function_avg(results[0].size(), 0.0);
		for (int m = 0; m < results.size(); ++m) {
			for (int i = 0; i < function_avg.size(); ++i) {
				function_avg[i] += results[m][i] / results.size();
			}
		}
		return function_avg;
	}



};