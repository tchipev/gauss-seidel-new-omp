/*
 * GaussSeidel2D.cpp
 *
 *  Created on: Feb 9, 2020
 *      Author: tchipevn
 */

#include "GaussSeidel2D.h"
#include <iostream>


void GaussSeidel2D::run(int whichSolver, int numIterations, int vtkOutput) {
	std::cout << "hi from " << whichSolver << std::endl;

	for(int y=1; y< _ny-1; ++y) {
		for(int x=1; x < _nx-1; ++x) {
			process9(x,y);
		}
	}
}
