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
#include <array>

#include "Slice1D.h"

class GaussSeidel2D {
public:
	/**
	 * @vtkOutput if a value <= 0 is passed, don't write files
	 */
	GaussSeidel2D(std::array<int, 2> N, std::array<int, 2> T, int vtkOutput) :
		_nx(N[0]), _ny(N[1]), _Tx(T[0]), _Ty(T[1]), _slice1dY(_Ty, _ny-2) {

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

	void printSchemeInfo(int solver) const;

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
	int _Tx, _Ty;

	std::vector<double> _values;

	Slice1D _slice1dY;

private:
	// **** traversals **** //

	// basic
	double basicTraversal() {
		double sumDiff2 = 0.0;
		#pragma omp parallel for reduction(+:sumDiff2) collapse(2)
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
		#pragma omp parallel for reduction(+:sumDiff2) collapse(2)
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

	double c04_hcpTraversal2() {
		double sumDiff2 = 0.0;

		int lookup[3] = {3, -3, -3};

		#pragma omp parallel reduction(+:sumDiff2)
		for (int col = 0; col < 3; ++col) {

			# pragma omp for collapse(3)
			for (int yOut = 1; yOut < _ny-1; yOut+=2) {
				for (int yIn = 0; yIn < 2; ++yIn) {
					for (int x = col * 2 + 1; x < _nx - 1 + 3; x += 6) {

						int yAccess = yIn + yOut;

						if (yAccess >= _ny-1)
							continue;

						int xAccess = x + yIn * lookup[col];

						if (xAccess >= _nx-1)
							continue;

						if(xAccess > 0)
							sumDiff2 += process9_residual(xAccess,yAccess);
						if(xAccess + 1 < _nx-1)
							sumDiff2 += process9_residual(xAccess+1,yAccess);
					}
				}
			}

			if(col < 2) {
				#pragma omp barrier
			}
		}
		return sumDiff2;
	}

	double c04_hcpTraversal3() {
		double sumDiff2 = 0.0;

		#pragma omp parallel reduction(+:sumDiff2)
		for (int col = 0; col < 3; ++col) {

			#pragma omp for collapse(2)
			for (int y = 1; y < _ny-1; ++y) {
				for (int x = col * 2 + 1; x < _nx-1 + 3; x += 6) {
					int xAccess = x + ((y-1)%2) * -3;

					if (xAccess > 0 and xAccess < _nx-1)
						sumDiff2 += process9_residual(xAccess,y);

					xAccess++;
					if (xAccess > 0 and xAccess < _nx-1)
						sumDiff2 += process9_residual(xAccess,y);
				}
			}
			if(col < 2) {
				#pragma omp barrier
			}
		}

		return sumDiff2;
	}

	double sli_along_y() {

		double sumDiff2 = 0.0;

		// determine max num threads that can be used
		int actualThreads = _slice1dY.getActualThreads();

		#pragma omp parallel num_threads(actualThreads) reduction(+:sumDiff2)
		{
			int myId = omp_get_thread_num();

			int myStartY = _slice1dY.getStart(myId);
			int myEndY = _slice1dY.getEnd(myId);

			_slice1dY.acquireLock(myId, Slice1D::MY_LOCK);

			// process front boundary
			int y = myStartY;
			for (int x = 1; x < _nx -1 ; ++x) {
				sumDiff2 += process9_residual(x,y+1);
			}

			_slice1dY.releaseLock(myId, Slice1D::MY_LOCK);

			// process bulk
			for (int y = myStartY + 1; y < myEndY-1; ++y) {
				for (int x = 1; x < _nx-1; ++x) {
					sumDiff2 += process9_residual(x,y+1);
				}
			}

			// acquire lock
			_slice1dY.acquireLock(myId, Slice1D::NEXT_LOCK);

			// process back boundary
			y = myEndY-1;
			for (int x = 1; x < _nx -1 ; ++x) {
				sumDiff2 += process9_residual(x,y+1);
			}

			// release lock
			_slice1dY.releaseLock(myId, Slice1D::NEXT_LOCK);

		}

		return sumDiff2;
	}

	double xFront(int xStart, int xEnd, int y) {
		double sumDiff2 = 0.0;
		for (int x = xStart; x < xEnd; ++x) {
			sumDiff2 += process9_residual(x, y);
		}
		return sumDiff2;
	}

	double yFront(int yStart, int yEnd, int x) {
		double sumDiff2 = 0.0;
		for (int y = yStart; y < yEnd; ++y) {
			sumDiff2 += process9_residual(x, y);
		}
		return sumDiff2;
	}

	double xyBlock(int xStart, int xEnd, int yStart, int yEnd) {
		double sumDiff2 = 0.0;
		for (int y = yStart; y < yEnd; ++y) {
			for (int x = xStart; x < xEnd; ++x) {
				sumDiff2 += process9_residual(x, y);
			}
		}
		return sumDiff2;
	}

	double sli_blk() {
		double sumDiff2 = 0.0;

		#pragma omp parallel reduction(+:sumDiff2)
		{
			int xStart, xEnd, yStart, yEnd;

			// acquire myH, myV

			// blue front
			sumDiff2 += xFront(xStart, xEnd/2, yStart);

			// release myH

			// blue block
			sumDiff2 += xyBlock(xStart+1, xEnd, yStart, yEnd);

			// acquire nbH

			// green front

		}
		return sumDiff2;
	}

};

#endif /* GAUSSSEIDEL2D_H_ */
