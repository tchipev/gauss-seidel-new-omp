/*
 * GaussSeidel2D.h
 *
 *  Created on: Feb 9, 2020
 *      Author: tchipevn
 */

#ifndef GAUSSSEIDEL2D_H_
#define GAUSSSEIDEL2D_H_

#include <vector>

class GaussSeidel2D {
public:
	GaussSeidel2D(int nx, int ny) :_nx(nx), _ny(ny) {
		_values.reserve(_nx * _ny);
		_values.resize(_nx * _ny);
	}
	~GaussSeidel2D() {
	}

	void run(int whichSolver);

	int getNx() const {
		return _nx;
	}

	int getNy() const {
		return _ny;
	}

private:
	// methods
	int ind(int x, int y) const {
		return x + _nx * y;
	}

	int numValues() const {
		return _nx * _ny;
	}

private:
	// fields
	int _nx, _ny;

	std::vector<double> _values;

};

#endif /* GAUSSSEIDEL2D_H_ */
