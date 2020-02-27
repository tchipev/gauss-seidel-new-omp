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
#include <getopt.h>

using namespace std;

struct Parameters {
	int _dimension;

	int _solver;
	int _iterations;
	int _vtkOutput;

	/* grid points */
	int _Nx, _Ny, _Nz;

	/* threads per direction */
	int _Tx, _Ty, _Tz;
};

Parameters readOptions(int argc, char * argv[]);

int main(int argc, char* argv[]) {

	Parameters p;
	try {
		p = readOptions(argc, argv);
	} catch(std::exception& e){
		return EXIT_FAILURE;
	}

	assert(p._dimension == 2 or p._dimension == 3);

	if (p._vtkOutput <= 0) {
		cout << "Program will not write VTK output and will write std out at every 10% of the progress." << endl;
	} else {
		cout << "Program will write VTK output every " << p._vtkOutput << " iterations and std out for every iteration." << endl;
	}

	if (p._dimension == 2) {
		std::array<int, 2> N = {p._Nx, p._Ny};
		std::array<int, 2> T = {p._Tx, p._Ty};

		GaussSeidel2D gs2D = GaussSeidel2D(N, T, p._vtkOutput);

		gs2D.run(p._solver, p._iterations, p._vtkOutput);
	} else if (p._dimension == 3) {
		std::array<int, 3> N = {p._Nx, p._Ny, p._Nz};
		std::array<int, 3> T = {p._Tx, p._Ty, p._Tz};

		GaussSeidel3D gs3D = GaussSeidel3D(N, T, p._vtkOutput);
		gs3D.run(p._solver, p._iterations, p._vtkOutput);
	}

	cout << "Exiting with success." << endl;
	return EXIT_SUCCESS;
}

Parameters readOptions(int argc, char * argv[]) {
	Parameters p;

	int opt;

	while ((opt = getopt(argc, argv, "d:s:i:o:x:y:z:u:v:w:?")) != -1) {
		switch (opt) {
		case 'd':
			p._dimension = atoi(optarg);
			cout << "Dimension: " << p._dimension << endl;
			break;
		case 's':
			p._solver = atoi(optarg);
			cout << "Solver: " << p._solver << endl;
			break;
		case 'i':
			p._iterations = atoi(optarg);
			cout << "Iterations: " << p._iterations << endl;
			break;
		case 'o':
			p._vtkOutput = atoi(optarg);
			cout << "VTK output: " << p._vtkOutput << endl;
			break;
		case 'x':
			p._Nx = atoi(optarg);
			cout << "Nx: " << p._Nx << endl;
			break;
		case 'y':
			p._Ny = atoi(optarg);
			cout << "Ny: " << p._Ny << endl;
			break;
		case 'z':
			p._Nz = atoi(optarg);
			cout << "Nz: " << p._Nz << endl;
			break;
		case 'u':
			p._Tx = atoi(optarg);
			cout << "Threads x: " << p._Tx << endl;
			break;
		case 'v':
			p._Ty = atoi(optarg);
			cout << "Threads y: " << p._Ty << endl;
			break;
		case 'w':
			p._Tz = atoi(optarg);
			cout << "Threads z: " << p._Tz << endl;
			break;
		case ':':
			cerr << "Missing arguments" << endl;
			throw std::exception();
		case '?':
			cerr << "Need 10ish input arguments to run." << endl;
			throw std::exception();
		}
	}

	return p;
}
