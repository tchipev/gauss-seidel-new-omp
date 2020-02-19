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
	/**
	 * @vtkOutput if a value <= 0 is passed, don't write files
	 */
	GaussSeidel2D(int nx, int ny, int vtkOutput) :_nx(nx), _ny(ny) {
		_values.reserve(_nx * _ny);
		_values.resize(_nx * _ny, 0.0);
		boundaryConditions(vtkOutput);
		initialConditions(vtkOutput);
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
	void boundaryConditions(int vtkOutput);

	void initialConditions(int vtkOutput);

	void printInfo(int iterations, double timeSec) const;

	void writeVTK(int iteration);

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

private:
	// **** traversals **** //

	// basic
	double basicTraversal() {
		double sumDiff2 = 0.0;
		for(int y=1; y < _ny-1; ++y) {
			for(int x=1; x < _nx-1; ++x) {
				// with residuum calculation
				double vold = val(x,y);

				process9(x,y);

				double diff = val(x,y) - vold;
				double diff2 = diff * diff;

				sumDiff2 += diff2;
			}
		}
		return sumDiff2;
	}

	// slow
	double slowTraversal() {
		double sumDiff2 = 0.0;

		for(int x=1; x < _nx-1; ++x) {
			for(int y=1; y < _ny-1; ++y) {
				// with residuum calculation
				double vold = val(x,y);

				process9(x,y);

				double diff = val(x,y) - vold;
				double diff2 = diff * diff;

				sumDiff2 += diff2;
			}
		}
		return sumDiff2;
	}

};

#endif /* GAUSSSEIDEL2D_H_ */
