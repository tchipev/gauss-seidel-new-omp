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
		_values.resize(_nx * _ny, 0.0);
		boundaryConditions();
		initialConditions();
	}
	~GaussSeidel2D() {
	}

	void run(int whichSolver, int numIterations, int vtkOutput);

	int getNx() const {
		return _nx;
	}

	int getNy() const {
		return _ny;
	}

private:
	void boundaryConditions();
	void initialConditions();

	void writeVTK(int iteration);

	// methods
	int ind(int x, int y) const {
		return x + _nx * y;
	}

	int numValues() const {
		return _nx * _ny;
	}

	void process9(int x, int y) {
		double diagSum = val(x-1,y-1)
					   + val(x+1,y-1)
					   + val(x-1,y+1)
					   + val(x+1,y+1);
		double lineSum = val(x,y-1)
					   + val(x,y+1)
					   + val(x-1,y)
					   + val(x+1,y);
		val(x,y) = ( (diagSum) + 4.0 * (lineSum) ) * 0.05;
	}

	double & val(int x, int y) {
		return _values[ind(x,y)];
	}

private:
	// fields
	int _nx, _ny;

	std::vector<double> _values;

};

#endif /* GAUSSSEIDEL2D_H_ */
