//============================================================================
// Name        : another-gauss-seidel.cpp
// Author      : Nikola Tchipev
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
//============================================================================

#include "GaussSeidel2D.h"
#include "GaussSeidel3D.h"
#include <stdlib.h>
#include <cstdlib>
#include <iostream>
#include <cassert>

using namespace std;

int main(int argc, char* argv[]) {

	if (argc < 2) {
		cout << "Expecting 6 arguments for 2D, 7 arguments for 3D." << endl;
		return EXIT_FAILURE;
	}

	int dimension=atoi(argv[1]);
	assert(dimension == 2 or dimension == 3);
	cout << "Dimension: " << dimension << endl;

	if ((dimension == 2 and argc != 7)
			or (dimension == 3 and argc != 8)) {
		cout << "Wrong number of arguments. Expecting 6 in 2D, 7 in 3D, found: " << argc - 1 << endl;
		cout << "exiting." << endl;
		return EXIT_FAILURE;
	}


	int whichSolver = atoi(argv[2]);
	cout << "using solver: " << whichSolver << endl;
	int numIterations = atoi(argv[3]);
	cout << "number of iterations: " << numIterations << endl;
	int vtkOutput = atoi(argv[4]);
	cout << "vtk output parameter: " << vtkOutput << endl;

	int nx = atoi(argv[5]);
	cout << "nx: " << nx << endl;
	int ny = atoi(argv[6]);
	cout << "ny: " << ny << endl;
	int nz;
	if (dimension == 3) {
		nz = atoi(argv[7]);
		cout << "nz: " << nz << endl;
	}

	if (vtkOutput <= 0) {
		cout << "Program will not write VTK output and will write std out at every 10% of the progress." << endl;
	} else {
		cout << "Program will write VTK output every " << vtkOutput << " iterations and std out for every iteration." << endl;
	}

	if (dimension == 2) {
		cout << "Running in 2D." << endl;
		GaussSeidel2D gs2D = GaussSeidel2D(nx, ny, vtkOutput);

		gs2D.run(whichSolver, numIterations, vtkOutput);
	} else if (dimension == 3) {
		cout << "Running in 3D." << endl;
		GaussSeidel3D gs3D = GaussSeidel3D(nx, nz, nz, vtkOutput);
		gs3D.run(whichSolver, numIterations, vtkOutput);
	}

	cout << "Exiting with success." << endl;
	return EXIT_SUCCESS;
}
