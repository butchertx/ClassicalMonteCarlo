#pragma once
#include <assert.h>

#include <tuple>
#include <vector>
#include <complex>

#include "cmctype.h"
#include "Lattice.h"
#include "RandomEngine.h"
#include "fftw_helper.h"


const double TAU = 2.0 * 3.14159265358979323846;

struct spin1 {

public:

	//holds the value (+-1) and the position indices of a given spin
	int s;
	int site;
};

class spin2 {

public:

	//holds the x- and y- components of a given spin
	cmctype::vec2<double> s;
	int site;

	spin2() {
		s.x = 1.0;
		s.y = 0.0;
		site = 0;
	}

	spin2(cmctype::vec2<double> s_) : s(s_), site(0) {};

	spin2(cmctype::vec2<double> s_, int site_) : s(s_), site(site_) {};

	spin2(double sx_, double sy_) : s(sx_, sy_), site(0) {};

	spin2(double sx_, double sy_, int site_) : s(sx_, sy_), site(site_) {};

	double dot(cmctype::vec2<double> v_) {
		return s * v_;
	}

	void normalize() {
		double norm = sqrt(s.x * s.x + s.y * s.y);
		if (std::abs(norm) < 1e-10) {
			s.x = 1.0;
			s.y = 0.0;
		}
		else {
			s.x = s.x / norm;
			s.y = s.y / norm;
		}
	}

	spin2 operator+(const spin2& s_) const {
		return spin2(s_.s + this->s, this->site);
	}

	spin2 operator-(const spin2& s_) const {
		return spin2(this->s - s_.s, this->site);
	}

	spin2 operator*(const double& c) const {
		return spin2((this->s)*c, this->site);
	}

	double operator*(const spin2& s_) const {
		return this->s * s_.s;
	}

	void operator+=(const spin2& s_) {
		this->s = s_.s + this->s;
	}
};

class spin3 {

public:

	//holds the x-, y-, and z- components of a given spin
	cmctype::vec3<double> s;
	int site;

	spin3() {
		s.x = 1.0;
		s.y = 0.0;
		s.z = 0.0;
		site = 0;
	}

	spin3(cmctype::vec3<double> s_) : s(s_), site(0) {};

	spin3(cmctype::vec3<double> s_, int site_) : s(s_), site(site_) {};

	spin3(double sx_, double sy_, double sz_) : s(sx_, sy_, sz_), site(0) {};

	spin3(double sx_, double sy_, double sz_, int site_) : s(sx_, sy_, sz_), site(site_) {};

	void normalize() {
		double norm = sqrt(s.x * s.x + s.y * s.y + s.z * s.z);
		if (std::abs(norm) < 1e-10) {
			s.x = 1.0;
			s.y = 0.0;
			s.z = 0.0;
		}
		else {
			s.x = s.x / norm;
			s.y = s.y / norm;
			s.z = s.z / norm;
		}
	}

	spin3 operator+(const spin3& s_) const {
		return spin3(s_.s + this->s, this->site);
	}

	spin3 operator*(const double& c) const {
		return spin3((this->s) * c, this->site);
	}

	double operator*(const spin3& s_) const {
		return this->s * s_.s;
	}

	void operator+=(const spin3& s_) {
		this->s = s_.s + this->s;
	}
};

template<class spintype>
class NRotorLattice {

protected:

	Lattice& lattice;
	const int dim = 3, N;
	const vec3<int> L;
	std::vector<spintype> spins;

	virtual void get_std_vector(std::vector<std::vector<double>>& spin_component_vectors) = 0;

public:

	NRotorLattice(Lattice& lattice_) : lattice(lattice_), L(lattice_.get_L()), N(lattice_.get_N()) {
		spins = std::vector<spintype>(N);
		for (int i = 0; i < N; ++i) {
			spins[i].site = i;
		}
	}

	int size() {	
		return spins.size();
	}

	void set(int site, spintype newspin) {
		spins[site] = newspin;
		spins[site].site = site;
	}

	spintype get(int site) {
		return spins[site];
	}

	virtual void randomize(RandomEngine* rand_p) = 0;

	virtual double calc_mag2() = 0;

	virtual double calc_chi_q(vec3<double> Q) = 0;

	virtual std::vector<double> calc_correlation() = 0;

};

class IsingLattice {
	//n = 1 rotor lattice in arbitrary dimension
	int dim;
	std::vector<int> L;
	std::vector<spin1> spins;
};

class XYLattice: public NRotorLattice<spin2> {
	//n = 2 rotor lattice in arbitrary dimension

	void get_std_vector(std::vector<std::vector<double>>& spin_component_vectors) {
		assert(spin_component_vectors.size() == 2);
		if (spin_component_vectors[0].size() != N || spin_component_vectors[1].size() != N) {
			spin_component_vectors[0].resize(N);
			spin_component_vectors[1].resize(N);
		}
		for (int i = 0; i < N; ++i) {
			spin_component_vectors[0][i] = spins[i].s.x;
			spin_component_vectors[1][i] = spins[i].s.y;
		}
	}

public:

	XYLattice(Lattice& lattice_) : NRotorLattice<spin2>(lattice_) {};

	void randomize(RandomEngine* rand_p) {
		double angle;
		spin2 R;
		for (int site = 0; site < spins.size(); ++site) {
			double angle = rand_p->get_rand_prob() * TAU;
			R = spin2(cos(angle), sin(angle));
			R.site = site;
			spins[site] = R;
		}
	}

	double calc_mag2() {
		vec2<double> result(0.0, 0.0);
		for (int i = 0; i < spins.size(); ++i) {
			result.x += spins[i].s.x;
			result.y += spins[i].s.y;
		}
		return (result * result) / ((double)N * N);
	}

	double calc_chi_q(vec3<double> Q) {
		vec2<std::complex<double>> result({ 0.0, 0.0 }, { 0.0, 0.0 });
		double phi = 0;
		for (int i = 0; i < spins.size(); ++i) {
			phi = lattice.get_coordinate(i) * Q;
			result.x += spins[i].s.x * std::complex<double>(cos(phi), sin(phi));
			result.y += spins[i].s.y * std::complex<double>(cos(phi), sin(phi));
		}
		return (std::abs(result.x) * std::abs(result.x) + std::abs(result.y) * std::abs(result.y)) / ((double)N * N);
	}

	std::vector<double> calc_correlation() {
		std::vector<std::vector<double>> sxy = { std::vector<double>(N, 0.0), std::vector<double>(N, 0.0) };
		get_std_vector(sxy);

		std::vector<double> corr_sx = calc_correlation_no_wisdom(sxy[0], L.x, L.y, L.z);
		std::vector<double> corr_sy = calc_correlation_no_wisdom(sxy[1], L.x, L.y, L.z);

		std::vector<double> result(N, 0.0);
		for (int i = 0; i < N; ++i) {
			result[i] += (corr_sx[i] + corr_sy[i]) / N;
		}
		return result;
	}

};

class XYZLattice : public NRotorLattice<spin3> {
	//n = 3 rotor lattice in arbitrary dimension

	void get_std_vector(std::vector<std::vector<double>>& spin_component_vectors) {
		assert(spin_component_vectors.size() == 3);
		if (spin_component_vectors[0].size() != N || spin_component_vectors[1].size() != N || spin_component_vectors[2].size() != N) {
			spin_component_vectors[0].resize(N);
			spin_component_vectors[1].resize(N);
			spin_component_vectors[2].resize(N);
		}
		for (int i = 0; i < N; ++i) {
			spin_component_vectors[0][i] = spins[i].s.x;
			spin_component_vectors[1][i] = spins[i].s.y;
			spin_component_vectors[2][i] = spins[i].s.z;
		}
	}

public:

	XYZLattice(Lattice& lattice_) : NRotorLattice<spin3>(lattice_) {};

	void randomize(RandomEngine* rand_p) {
		double angle;
		spin3 R;
		for (int site = 0; site < spins.size(); ++site) {
			double phi = rand_p->get_rand_prob() * TAU;
			double theta = rand_p->get_rand_prob() * 0.5 * TAU;
			R = spin3(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
			R.site = site;
			spins[site] = R;
		}
	}

	double calc_mag2() {
		vec3<double> result(0.0, 0.0, 0.0);
		for (int i = 0; i < spins.size(); ++i) {
			result.x += spins[i].s.x;
			result.y += spins[i].s.y;
			result.z += spins[i].s.z;
		}
		return (result * result) / ((double)N * N);
	}

	double calc_chi_q(vec3<double> Q) {
		vec3<std::complex<double>> result({ 0.0, 0.0 }, { 0.0, 0.0 }, { 0.0, 0.0 });
		double phi = 0;
		for (int i = 0; i < spins.size(); ++i) {
			phi = lattice.get_coordinate(i) * Q;
			result.x += spins[i].s.x * std::complex<double>(cos(phi), sin(phi));
			result.y += spins[i].s.y * std::complex<double>(cos(phi), sin(phi));
			result.z += spins[i].s.z * std::complex<double>(cos(phi), sin(phi));
		}
		return (std::abs(result.x) * std::abs(result.x) + std::abs(result.y) * std::abs(result.y) + std::abs(result.z) * std::abs(result.z)) / ((double)N * N);
	}

	std::vector<double> calc_correlation() {
		std::vector<std::vector<double>> sxyz = { std::vector<double>(N, 0.0), std::vector<double>(N, 0.0), std::vector<double>(N, 0.0) };
		get_std_vector(sxyz);

		std::vector<double> corr_sx = calc_correlation_no_wisdom(sxyz[0], L.x, L.y, L.z);
		std::vector<double> corr_sy = calc_correlation_no_wisdom(sxyz[1], L.x, L.y, L.z);
		std::vector<double> corr_sz = calc_correlation_no_wisdom(sxyz[2], L.x, L.y, L.z);

		std::vector<double> result(N, 0.0);
		for (int i = 0; i < N; ++i) {
			result[i] += (corr_sx[i] + corr_sy[i] + corr_sz[i]) / N;
		}
		return result;
	}

};