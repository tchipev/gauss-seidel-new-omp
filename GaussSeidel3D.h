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
	// UNSAFE OMP parallelization!
	double basicTraversal() {
		double sumDiff2 = 0.0;
		#pragma omp parallel for reduction(+:sumDiff2) collapse(3)
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
	// UNSAFE OMP parallelization!
	double slowTraversal() {
		double sumDiff2 = 0.0;
		#pragma omp parallel for reduction(+:sumDiff2) collapse(3)
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
		double sumDiff2 = 0.0;

		#pragma omp parallel reduction(+:sumDiff2)
		for (int colour = 0; colour < 8; ++colour) {
			int startX = colour % 2 + 1;
			int startY = (colour / 2) % 2 + 1;
			int startZ = colour / 4 + 1;

			#pragma omp for nowait collapse(3)
			for (int z=startZ; z < _nz-1; z+=2) {
				for (int y=startY; y < _ny-1; y+=2) {
					for (int x=startX; x < _nx-1; x+=2) {
						sumDiff2 += process27_residual(x,y,z);
					}
				}
			}

			if (colour < 7) {
				#pragma omp barrier
			}
		} /* end of for, end of parallel */

		return sumDiff2;
	}

	int myMod1(int x) {
		int mod = x % 12;
		int ret = mod;
		switch(mod) {
		case 10: /* FALLTHRU */
		case 11:
			ret = mod - 12;
			break;
		default:
			break;
		}
		return ret;
	}

	int myMod2(int x) {
		int mod = x % 12;
		return mod < 10 ? mod : mod - 12;
	}

	int myMod3(int x) {
	int mod = x % 12;
	int div = mod / 10;
	return mod - div * 12;
}

#if 0
	// c04 data:
	std::vector<std::pair<int,int>> _c04_indices;

	void prepareIndices() {
		for (int col = 0; col < 4; ++col) {

			// TODO: account for boundaries
			for (int z = 1; z < _nz-1; ++z) {

				if (z % 2 == 0) {
					// displaced plane
				}
				int startY;

				for (int y = 1; y < _ny-1; ++y) {
					int startX;
					if (y % 2 == 0) {
						startX = col * 2;
					} else {
						startX =
					}

					for (int x = startX; x < _nx; ++x) {

					}
				}

			}

		}
	}


	// c04_hcp
	double c04_hcpTraversal() {
		double sumDiff2 = 0.0;

		for (int colour = 0; colour < 4; ++colour) {

	//		int startX[2][4][3] = {{{0,4,8},{3,7,-1},{6,-2,2}, {9,1,5}}};
			for (int z = 1; z < _nz-1 ; ++z) {
				for (int y = 1; y < _ny-1; y+=2) {

					bool secondY = y+1 < _ny-1;

					int startX = ((y%6)/2)*4 + (z % 2)*4 + colour * 3;
					startX = myMod2(startX);

					for (int x = startX; x < _nx-1; x+=12) {

					}
				}
			}
		}

		return sumDiff2;
	}
#else
	double c04_hcpTraversal() {
		double sumDiff2 = 0.0;

//		#pragma omp parallel reduction(+:sumDiff2)
		for (int col = 0; col < 4; ++col) {

//			#pragma omp for collapse(3)
			for (int z = 1; z < _nz-1; ++z) {
				for (int y = 1; y < _ny-1 + 3; y += 2) {
					for (int x = col * 3 + 1; x < _nx-1 + 8; x += 12) {

						int yAccess = y + ((z-1)%2) * -3;

						int xAccess = x + (((y-1)/2)%3) * -4;

						for(int j = 0; j < 2; ++j) {
							int Y = yAccess + j;
							if (Y > 0 and Y < _ny-1) {
								for (int i = 0; i < 3; ++i) {
									int X = xAccess + i;
									if (X > 0 and X <_nx-1) {
										sumDiff2 += process27_residual(X, Y, z);
									}
								}
							} // end i loop
						} // end j loop
					} // end x loop
				} // end y loop
			} // end z loop
			if(col < 3) {
//				#pragma omp barrier
			}

			writeVTK(1000 + col);
		} // end col
		return sumDiff2;
	}
#endif

};

#endif /* GAUSSSEIDEL3D_H_ */
