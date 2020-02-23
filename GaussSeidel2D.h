/*
 * GaussSeidel2D.h
 *
 *  Created on: Feb 9, 2020
 *      Author: tchipevn
 */

#ifndef GAUSSSEIDEL2D_H_
#define GAUSSSEIDEL2D_H_

#include <vector>
#include <omp.h>

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

private:
	void boundaryConditions(int vtkOutput);

	void initialConditions(int vtkOutput);

	void printInfo(int iterations, double timeSec) const;

	void writeVTK(int iteration);

	int ind(int x, int y) const {
		return x + _nx * y;
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

	double process9_residual(int x, int y) {
		double vold = val(x,y);

		process9(x,y);

		double diff = val(x,y) - vold;
		double diff2 = diff * diff;

		return diff2;
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
				sumDiff2 += process9_residual(x,y);
			}
		}
		return sumDiff2;
	}

	// slow
	double slowTraversal() {
		double sumDiff2 = 0.0;

		for(int x=1; x < _nx-1; ++x) {
			for(int y=1; y < _ny-1; ++y) {
				sumDiff2 += process9_residual(x,y);
			}
		}
		return sumDiff2;
	}

	// c08
	double c08Traversal() {
		double sumDiff2 = 0.0;

		#pragma omp parallel reduction(+:sumDiff2)
		for (int colour = 0; colour < 4; ++colour) {

			int startX = colour % 2 + 1;
			int startY = colour / 2 + 1;

			#pragma omp for nowait collapse(2)
			for(int y=startY; y < _ny-1; y += 2) {
				for(int x=startX; x < _nx-1; x += 2) {
					sumDiff2 += process9_residual(x,y);
				}
			}
			if(colour < 3) {
				#pragma omp barrier
			}
		} /* end of for, end of parallel */

		return sumDiff2;
	}

	// c04_hcp
	double c04_hcpTraversal() {
		double sumDiff2 = 0.0;
		int startX[3][2] = {{4,1},{0,3},{2,5}};

		#pragma omp parallel reduction(+:sumDiff2)
		for (int colour = 0; colour < 3; ++colour) {

			#pragma omp for nowait
			for(int y=1; y < _ny-1; ++y) {
				for(int x=startX[colour][y%2]; x < _nx-1; x += 6) {
					if(x > 0)
						sumDiff2 += process9_residual(x,y);
					if(x + 1 < _nx-1)
						sumDiff2 += process9_residual(x+1,y);
				}
			}

			if(colour < 2) {
				#pragma omp barrier
			}
		} /* end of for, end of parallel */

		return sumDiff2;
	}

};

#endif /* GAUSSSEIDEL2D_H_ */
