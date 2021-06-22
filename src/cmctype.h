#pragma once
#include <vector>
#include <complex>
#include <iostream>
#include <string>
#include <fstream>
#include "json.hpp"

#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif

// for convenience
using json = nlohmann::json;

namespace cmctype {

	struct LatticeParams {

		std::string spintype, latticetype;
		std::vector<int> L;
		std::string pbc_type;

	};

	void to_json(json& j, const LatticeParams& p);

	void from_json(const json& j, LatticeParams& p);

	struct Interaction {

		std::string name;
		double strength;

	};

	void to_json(json& j, const Interaction& p);

	void from_json(const json& j, Interaction& p);

	struct ModelParams {

		std::string model_name;
		std::vector<Interaction> interactions;
		double beta;
		double T;

		double get_interaction(std::string interaction_name) {
			for (auto i : interactions) {
				if (i.name.compare(interaction_name) == 0) {
					return i.strength;
				}
			}
			std::cout << "No valid interaction found for interaction name: " << interaction_name << "\n";
		}

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

		void print();

		void set_parallel_params(int id);

	};

	MCParams read_params(json j);

	MCParams read_params(std::string infile_name);

	/*! \brief Small value for double comparisons
		*/
	const double EPSILON = 1e-12;

	/*! \brief 2-D Vector
	*/
	template <class T>
	class vec2 {
	public:
		T x, y;

		vec2() {}

		vec2(T x_in, T y_in) : x(x_in), y(y_in) {}

		vec2<T> operator+ (const vec2<T>& v) const {
			return vec2<T>(this->x + v.x, this->y + v.y);
		}

		vec2<T> operator- (const vec2<T>& v) const {
			return vec2<T>(this->x - v.x, this->y - v.y);
		}

		vec2<T> operator* (const double& c) const {
			return vec2<T>(this->x*c, this->y*c);
		}

		double operator* (const vec2<T>& v) const {
			return this->x * v.x + this->y * v.y;
		}

	};

	/*! \brief 3-D Vector with operations
	*/
	template <class T>
	class vec3 {
	public:
		T x, y, z;

		vec3() {}

		vec3(std::vector<T> v_in) {
			assert(v_in.size() == 3);
			x = v_in[0];
			y = v_in[1];
			z = v_in[2];
		}

		vec3(T x_in, T y_in, T z_in) : x(x_in), y(y_in), z(z_in) {}

		double angle_xy() {
			//get the angle in rad in the xy plane (from the x-axis)
			if (y == 0) {
				return x > 0 ? 0.0 : M_PI;
			}
			else {
				return atan2(y, x) < 0 ? atan2(y, x) + 2.0*M_PI : atan2(y, x);//[0,2pi]
			}
		}

		vec3 operator+ (const vec3<T>& v) const {
			return vec3(this->x + v.x, this->y + v.y, this->z + v.z);
		}

		vec3 operator* (const int& a) const {
			return vec3(this->x * a, this->y * a, this->z * a);
		}

		vec3 operator* (const double& a) const {
			return vec3(this->x * a, this->y * a, this->z * a);
		}

		double operator* (const vec2<T>& v) const {
			return this->x * v.x + this->y * v.y + this->z * v.z;
		}

		template <class R>
		R operator* (const vec3<R>& v) const {
			return (v.x * this->x + v.y * this->y + v.z * this->z);
		}

		vec3 operator% (const vec3<T>& v) const {
			//cross product
			return vec3(this->y * v.z - this->z * v.y, this->z * v.x - this->x * v.z, this->x * v.y - this->y * v.x);
		}

		vec3 operator- (const vec3<T>& v) const {
			return vec3(this->x - v.x, this->y - v.y, this->z - v.z);
		}

		double operator^ (const vec3<T>& v) const {
			//distance between the two vectors
			vec3 diff(this->x - v.x, this->y - v.y, this->z - v.z);
			return sqrt(diff * diff);
		}

		bool operator== (const vec3<T>& v) const {
			//compare vectors
			return (std::abs(this->x - v.x) < EPSILON
				&& std::abs(this->y - v.y) < EPSILON
				&& std::abs(this->z - v.z) < EPSILON);
		}

	};

	template <class T>
	std::string vec3str(vec3<T> vec_in) {
		std::stringstream ss;
		ss << "<" << vec_in.x << "," << vec_in.y << "," << vec_in.z << ">";
		return ss.str();
	}

	template <class T>
	void to_json(json& j, const vec3<T>& p) {
		j = json{ vec3str(p) };
	}

	template <class T>
	void from_json(const json& j, vec3<T>& p) {
		std::vector<T> vec(3);
		j.get_to<std::vector<T>>(vec);
		p.x = vec[0];
		p.y = vec[1];
		p.z = vec[2];
	}

}