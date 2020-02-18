//============================================================================
// Name        : another-gauss-seidel.cpp
// Author      : Nikola Tchipev
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
//============================================================================

#include "GaussSeidel2D.h"
#include <stdio.h>
#include <stdlib.h>


int main(void) {
	puts("Hodor!");

	int nx = 10, ny = 10;
	int whichSolver = 0;
	int numIterations = 10;
	int vtkOutput = 1;

	GaussSeidel2D gs2D = GaussSeidel2D(nx, ny);

	gs2D.run(whichSolver, numIterations, vtkOutput);


	return EXIT_SUCCESS;
}
