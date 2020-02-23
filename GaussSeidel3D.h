/*
 * GaussSeidel3D.h
 *
 *  Created on: Feb 23, 2020
 *      Author: tchipevn
 */

#ifndef GAUSSSEIDEL3D_H_
#define GAUSSSEIDEL3D_H_

#include <vector>
#include <omp.h>

class GaussSeidel3D {
public:
	GaussSeidel3D(int nx, int ny, int nz, int vtkOutput) :_nx(nx), _ny(ny), _nz(nz) {
		int numValues = _nx * _ny * _nz;
		_values.reserve(numValues);
		_values.resize(numValues, 0.0);
		boundaryConditions(vtkOutput);
		initialConditions(vtkOutput);
	}
	~GaussSeidel3D() {
		// TODO Auto-generated destructor stub
	}

	void run(int whichSolver, int numIterations, int vtkOutput);

private:
	// methods
	void boundaryConditions(int vtkOutput);

	void initialConditions(int vtkOutput);

	void printInfo(int iterations, double timeSec) const;

	void writeVTK(int iteration);

	int ind(int x, int y, int z) const {
		return x + _nx * y + _nx * _ny * z;
	}

	void process27(int x, int y, int z) {
		double faceSum = val(x,y,z-1)
					   + val(x,y,z+1)
					   + val(x,y-1,z)
					   + val(x,y+1,z)
					   + val(x-1,y,z)
					   + val(x+1,y,z);

		double edgeSum = val(x,y-1,z-1)
					   + val(x,y+1,z-1)
					   + val(x-1,y,z-1)
					   + val(x+1,y,z-1)

					   + val(x-1,y-1,z)
					   + val(x-1,y+1,z)
					   + val(x+1,y-1,z)
					   + val(x+1,y+1,z)

					   + val(x,y-1,z+1)
					   + val(x,y+1,z+1)
					   + val(x-1,y,z+1)
					   + val(x+1,y,z+1);

		double cornerSum = val(x-1,y-1,z-1)
						 + val(x-1,y+1,z-1)
						 + val(x-1,y+1,z-1)
						 + val(x+1,y+1,z-1)

						 + val(x-1,y-1,z+1)
						 + val(x-1,y+1,z+1)
						 + val(x-1,y+1,z+1)
						 + val(x+1,y+1,z+1);
		val(x,y,z) = (cornerSum + 3.0 * edgeSum + 14.0 * faceSum) * 0.0078125;
	}

	double process27_residual(int x, int y, int z) {
		double vold = val(x,y,z);

		process27(x,y,z);

		double diff = val(x,y,z) - vold;
		double diff2 = diff * diff;

		return diff2;
	}

	double & val(int x, int y, int z) {
		return _values[ind(x,y,z)];
	}

private:
	// fields
	int _nx, _ny, _nz;

	std::vector<double> _values;

private:
	// **** traversals **** //

	// basic
	double basicTraversal() {
		double sumDiff2 = 0.0;
		for (int z = 1; z < _nz-1; ++z) {
			for (int y = 1; y < _ny-1; ++y) {
				for (int x = 1; x < _nx-1; ++x) {
					sumDiff2 += process27_residual(x,y,z);
				}
			}
		}
		return sumDiff2;
	}

	// slow
	double slowTraversal() {
		double sumDiff2 = 0.0;
		for (int x = 1; x < _nx-1; ++x) {
			for (int y = 1; y < _ny-1; ++y) {
				for (int z = 1; z < _nz-1; ++z) {
					sumDiff2 += process27_residual(x,y,z);
				}
			}
		}
		return sumDiff2;
	}

	// c08
	double c08Traversal() {
		return 0.0;
	}

	// c04_hcp
	double c04_hcpTraversal() {
		return 0.0;
	}
};

#endif /* GAUSSSEIDEL3D_H_ */
