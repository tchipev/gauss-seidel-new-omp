//============================================================================
// Name        : another-gauss-seidel.cpp
// Author      : Nikola Tchipev
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
//============================================================================

#include "GaussSeidel2D.h"
#include <stdlib.h>
#include <cstdlib>
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {

	if (argc != 6) {
		cout << "Wrong number of arguments. Expected 5, found: " << argc -1 << endl;
		cout << "exiting." << endl;
		return EXIT_FAILURE;
	}

	int nx = atoi(argv[1]);
	int ny = atoi(argv[2]);
	int whichSolver = atoi(argv[3]);
	int numIterations = atoi(argv[4]);
	int vtkOutput = atoi(argv[5]);

	GaussSeidel2D gs2D = GaussSeidel2D(nx, ny);

	gs2D.run(whichSolver, numIterations, vtkOutput);

	cout << "Exiting with success." << endl;
	return EXIT_SUCCESS;
}
