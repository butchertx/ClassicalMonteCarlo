#include "MCResults.h"
#include <sys/stat.h> // stat
#include <errno.h>    // errno, ENOENT, EEXIST
#if defined(_WIN32)
#include <direct.h>   // _mkdir
#endif
#include <assert.h>
#include "fftw_helper.h"
#include <numeric>

bool isDirExist(const std::string& path)
{
#if defined(_WIN32)
	struct _stat info;
	if (_stat(path.c_str(), &info) != 0)
	{
		return false;
	}
	return (info.st_mode & _S_IFDIR) != 0;
#else 
	struct stat info;
	if (stat(path.c_str(), &info) != 0)
	{
		return false;
	}
	return (info.st_mode & S_IFDIR) != 0;
#endif
}

bool makePath(const std::string& path)
{
#if defined(_WIN32)
	int ret = _mkdir(path.c_str());
#else
	mode_t mode = 0755;
	int ret = mkdir(path.c_str(), mode);
#endif
	if (ret == 0)
		return true;

	switch (errno)
	{
	case ENOENT:
		// parent didn't exist, try to create it
	{
		int pos = path.find_last_of('/');
		if (pos == std::string::npos)
#if defined(_WIN32)
			pos = path.find_last_of('\\');
		if (pos == std::string::npos)
#endif
			return false;
		if (!makePath(path.substr(0, pos)))
			return false;
	}
	// now, try to create again
#if defined(_WIN32)
	return 0 == _mkdir(path.c_str());
#else 
	return 0 == mkdir(path.c_str(), mode);
#endif

	case EEXIST:
		// done!
		return isDirExist(path);

	default:
		return false;
	}
}

void MCResults::record(std::string obs, double measurement) {
	measurements[obs].push_back(measurement);
}

void MCResults::record(std::string obs, std::vector<double> measurement) {
	function_measurements[obs].push_back(measurement);

	if (!keep_functions &&
		function_measurements[obs].size() == function_bin_size) {
		bin_function(obs);
		function_measurements[obs].clear();
	}
}

void MCResults::bin_function(std::string obs) {
	assert(function_measurements[obs].size() == function_bin_size);

	//get measurements and create zero function to accumulate into
	std::vector<std::vector<double>> measures = function_measurements[obs];
	std::vector<double> cumul(measures[0].size(), 0.0);
	std::vector<double> cumul_sq(measures[0].size(), 0.0);

	for (int meas = 0; meas < measures.size(); ++meas) {
		for (int i = 0; i < measures[meas].size(); ++i) {
			cumul[i] += measures[meas][i] / (1.0 * function_bin_size);
			cumul_sq[i] += (measures[meas][i] * measures[meas][i]) / (1.0 * function_bin_size);
		}
	}

	function_bins[obs].push_back(cumul);
	function_sq_bins[obs].push_back(cumul_sq);
	function_counts[obs].push_back(function_bin_size);
}

std::vector<double> MCResults::get(std::string obs) {
	return measurements[obs];
}

std::vector<double> MCResults::get_function_average(std::string obs) {
	std::vector<std::vector<double>> results;
	if (keep_functions) {
		results = function_measurements[obs];
	}
	else {
		results = function_bins[obs];
	}
	std::vector<double> function_avg(results[0].size(), 0.0);
	for (int m = 0; m < results.size(); ++m) {
		for (int i = 0; i < function_avg.size(); ++i) {
			function_avg[i] += results[m][i] / results.size();
		}
	}
	return function_avg;
}

std::vector<double> MCResults::get_function_var(std::string obs) {
	std::vector<std::vector<double>> results;
	if (keep_functions) {
		results = function_measurements[obs];
	}
	else {
		results = function_sq_bins[obs];
	}
	std::vector<double> function_avg = get_function_average(obs);
	std::vector<double> function_var(results[0].size(), 0.0);

	for (int i = 0; i < function_var.size(); ++i) {
		for (int m = 0; m < results.size(); ++m) {
			function_var[i] += results[m][i] / results.size();
		}
		function_var[i] -= function_avg[i] * function_avg[i];
	}
	return function_var;
}

std::vector<double> MCResults::autocorrelation(std::vector<double> measurements_) {
	std::vector<double> result(measurements_.size() - 1, 0.0);
	std::vector<int> counts(measurements_.size() - 1, 0);

	double sum_m = std::accumulate(measurements_.begin(), measurements_.end()-1, 0.0);
	double sum_n = std::accumulate(measurements_.begin()+1, measurements_.end(), 0.0);
	double mean_m = sum_m / (measurements_.size() - 1);
	double mean_n = sum_n / (measurements_.size() - 1);

	for (int m = 0; m < measurements_.size()-1; ++m) {
		for (int n = m + 1; n < measurements_.size(); ++n) {
			result[n - m - 1] += measurements_[n] * measurements_[m];
			counts[n - m - 1] += 1;
		}
	}

	for (int i = 0; i < result.size(); ++i) {
		result[i] = result[i] / counts[i] - mean_m * mean_n;
	}

	return result;
}