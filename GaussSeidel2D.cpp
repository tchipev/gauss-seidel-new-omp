/*
 * GaussSeidel2D.cpp
 *
 *  Created on: Feb 9, 2020
 *      Author: tchipevn
 */

#include "GaussSeidel2D.h"
#include <iostream>
#include <fstream>   /* file io */
#include <stdlib.h>  /* srand */

using namespace std;

void GaussSeidel2D::run(int whichSolver, int numIterations, int vtkOutput) {
	std::cout << "hi from " << whichSolver << std::endl;

	for(int y=1; y< _ny-1; ++y) {
		for(int x=1; x < _nx-1; ++x) {
			process9(x,y);
		}
	}

	writeVTK(10);
}

void GaussSeidel2D::boundaryConditions() {
	// left boundary 0, right boundary 1
	for (int y = 0; y < _ny; ++y) {
		val(0,y) = 0.0;
		val(_nx-1,y) = 1.0;
	}
	// upper and lower boundaries linear interpolation
	double nx_inverse = 1.0 / (_nx-1);
	for (int x = 0; x < _nx; ++x) {
		val(x,0) = x * nx_inverse;
	}
	for (int x = 0; x < _nx; ++x) {
		val(x,_ny-1) = x * nx_inverse;
	}
	writeVTK(-2);
}

void GaussSeidel2D::initialConditions() {
	// pi + rand[0,1]
	// for fixed seed

	const int SEED = 42;
	srand(SEED);

	double rand_max_inv = 1.0 / RAND_MAX;

	for(int y=1; y< _ny-1; ++y) {
		for(int x=1; x < _nx-1; ++x) {
			val(x,y) = 3.14159265 + rand() * rand_max_inv;
		}
	}

	writeVTK(-1);
}

void GaussSeidel2D::writeVTK(int iteration) {
	ofstream myfile;
	string filename = "vtk/val_"+ to_string(iteration) + ".vtk";
	myfile.open (filename);
	myfile << "# vtk DataFile Version 2.0\n";
	myfile << "volume example " << endl;
	myfile << "ASCII " << endl;
	myfile << "DATASET STRUCTURED_POINTS " << endl;
	myfile << "DIMENSIONS " << _nx << " " << _ny << " 1" << "\n";
	myfile << "ASPECT_RATIO 1 1 1 " << endl;
	myfile << "ORIGIN 0 0 0 " << endl;
	myfile << "POINT_DATA " << _nx * _ny << "\n";
	myfile << "SCALARS volume_scalars double " << endl;
	myfile << "LOOKUP_TABLE default " << endl;

	for (int y = 0; y < _ny; ++y) {
		for (int x = 0; x < _nx; ++x) {
			myfile << val(x,y) << " ";
		}
		myfile << endl;
	}

	myfile.close();
}



