#pragma once
#include <vector>
#include <assert.h>
#include "MCParams.h"
#include "RandomEngine.h"
#include "fftw_helper.h"
#include <tuple>

const double TAU = 2.0 * 3.14159265358979323846;

template <unsigned int N = 2>
class Spin {

	double s[N];

};

struct spin1 {
	//holds the value (+-1) and the position indices of a given spin
	int s;
	int x;
	int y;
};

class spin2 {

public:

	//holds the x- and y- components of a given spin
	double sx, sy;
	int site;
	//int x;
	//int y;

	spin2() {
		sx = 1.0;
		sy = 0.0;
		site = 0;
	}

	spin2(double sx_, double sy_) : sx(sx_), sy(sy_), site(0) {};

	spin2(double sx_, double sy_, int site_) : sx(sx_), sy(sy_), site(site_) {};

	double dot(spin2 other) {
		return sx * other.sx + sy * other.sy;
	}

	void normalize() {
		double norm = sqrt(sx * sx + sy * sy);
		if (std::abs(norm) < 1e-10) {
			sx = 1.0;
			sy = 0.0;
		}
		else {
			sx = sx / norm;
			sy = sy / norm;
		}
	}

	spin2 operator+(const spin2& s) {
		return spin2(s.sx + this->sx, s.sy + this->sy);
	}

	spin2 operator*(const double& c) {
		return spin2(c*(this->sx), c*(this->sy));
	}

	void operator+=(const spin2& s) {
		this->sx += s.sx;
		this->sy += s.sy;
	}
};

struct spin3 {
	//holds the x- and y- and z- components of rotor and the position indices of a given spin
	double sx, sy, sz;
	int x;
	int y;
};

template<class spintype>
class NRotorLattice {

protected:

	int dim, N;
	std::vector<int> L;
	std::vector<spintype> spins;
	std::vector<std::vector<int>> rolled_indices;
	bool pbc = true;


	int unroll_index(std::vector<int> index) {
		assert(index.size() == L.size());
		int site = 0;
		int block = 1;
		for (int i = dim - 1; i >= 0; --i) {
			site += index[i] * block;
			block *= L[i];
		}
		return site;
	}

	std::vector<int> roll_index(int site) {
		return rolled_indices[site];
	}

public:

	NRotorLattice(int dim_, std::vector<int> L_) : dim(dim_), L(L_) {
		assert(dim == L.size());
		N = 1;
		for (int i = 0; i < dim; ++i) {
			N *= L[i];
		}
		spins = std::vector<spintype>(N);
		for (int i = 0; i < N; ++i) {
			spins[i].site = i;
		}
		if (dim == 3) {
			//unrolled spin list gives site = nx*Lz*Ly + ny*Lz + nz
			for (int site = 0; site < N; ++site) {
				rolled_indices.push_back({});
				rolled_indices[site].push_back((site / L[2]) / L[1]); //nx
				rolled_indices[site].push_back((site / L[2]) % L[1]); //ny
				rolled_indices[site].push_back(site % L[2]); //nz
			}
		}
	};

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

	spintype get(int site, std::vector<int> displacement) {
		
		assert(L.size() == displacement.size());
		if (pbc) {
			std::vector<int> rolled = rolled_indices[site];
			for (int i = 0; i < dim; ++i) {
				rolled[i] = (rolled[i] + displacement[i] + L[i]) % L[i];
			}
			return get(unroll_index(rolled));
		}
		
	}

	void get_std_vector(std::vector<double>& result) {
		if (result.size() != N) {
			result.resize(N);
		}
		for (int i = 0; i < N; ++i) {
			result[i] = (double)spins[i].s;
		}
	}

	void get_std_vector(std::vector<double>& result_x, std::vector<double>& result_y) {
		if (result_x.size() != N || result_y.size() != N) {
			result_x.resize(N);
			result_y.resize(N);
		}
		for (int i = 0; i < N; ++i) {
			result_x[i] = (double)spins[i].sx;
			result_y[i] = (double)spins[i].sy;
		}
	}

	virtual void randomize(RandomEngine* rand_p) {};

	virtual double calc_mag2() { return 0.0; }

};

class IsingLattice {
	//n = 1 rotor lattice in arbitrary dimension
	int dim;
	std::vector<int> L;
	std::vector<spin1> spins;
};

class XYLattice3D: public NRotorLattice<spin2> {
	//n = 2 rotor lattice in arbitrary dimension

public:

	XYLattice3D(std::vector<int> L_) : NRotorLattice<spin2>(3, L_) {
		assert(L.size() == 3);
		//spins = std::vector<spin2>(L[0] * L[1] * L[2]);//can go up to > 2000 x 1000 x 1000 without integer overflow problems
	};

	double calc_mag2() {
		spin2 result(0.0, 0.0);
		for (int i = 0; i < spins.size(); ++i) {
			result.sx += spins[i].sx;
			result.sy += spins[i].sy;
		}
		return ((double)result.sx * result.sx + result.sy * result.sy) / ((double) N * N);
	}

	double calc_chi_q(std::vector<double> Q) {
		assert(Q.size() == 3);
		spin2 result_r(0.0, 0.0), result_i(0.0, 0.0);
		double phi = 0;
		for (int i = 0; i < spins.size(); ++i) {
			//this is so hacky and bad and slow
			phi = roll_index(i)[0] * Q[0] + roll_index(i)[1] * Q[1] + roll_index(i)[2] * Q[2];
			result_r.sx += spins[i].sx * cos(phi);
			result_r.sy += spins[i].sy * cos(phi);
			result_i.sx += spins[i].sx * sin(phi);
			result_i.sy += spins[i].sy * sin(phi);
		}
		return ((double)result_r.sx * result_r.sx + result_r.sy*result_r.sy + result_i.sx*result_i.sx + result_i.sy*result_i.sy) / ((double)N * N);
	}

	std::vector<double> calc_correlation() {
		std::vector<double> sx(N, 0.0), sy(N, 0.0);
		get_std_vector(sx, sy);

		std::vector<double> corr_sx = calc_correlation_no_wisdom(sx, L[0], L[1], L[2]);
		std::vector<double> corr_sy = calc_correlation_no_wisdom(sy, L[0], L[1], L[2]);

		std::vector<double> result(N, 0.0);
		for (int i = 0; i < N; ++i) {
			result[i] += (corr_sx[i] + corr_sy[i]) / N;
		}
		return result;
	}

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

};

class XYZLattice {
	//n = 3 rotor lattice in arbitrary dimension
	int dim;
	std::vector<int> L;
	std::vector<spin3> spins;
};
