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
#include <chrono>    /* steady_clock */

using namespace std;
using namespace std::chrono;

void GaussSeidel2D::run(int whichSolver, int numIterations, int vtkOutput) {
	steady_clock::time_point t1 = steady_clock::now();

	int tenPercentOfIterations = numIterations / 10;
	if (tenPercentOfIterations == 0) {
		tenPercentOfIterations = 1;
	}

	for (int iteration = 1; iteration <= numIterations; ++iteration) {

		double sumDiff2;
		switch(whichSolver) {
		case 0:
			sumDiff2 = basicTraversal();
			break;
		case 1:
			sumDiff2 = slowTraversal();
			break;
		case 2:
			sumDiff2 = c08Traversal();
			break;
		case 3:
			sumDiff2 = c04_hcpTraversal();
			break;
		}

		if (vtkOutput > 0 or iteration % tenPercentOfIterations == 0) {
			cout << "iteration: " << iteration << " residuum square: " << sumDiff2 << endl;
		}

		if (vtkOutput > 0 and iteration % vtkOutput == 0) {
			cout << "Writing VTK output" << endl;
			writeVTK(iteration);
		}
	}

	if(vtkOutput > 0) {
		writeVTK(numIterations);
	}

	steady_clock::time_point t2 = steady_clock::now();

	// eclipse is not that smart?
	double time = duration_cast<duration<double>>(t2 - t1).count();

	printInfo(numIterations, time);
}

void GaussSeidel2D::boundaryConditions(int vtkOutput) {
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

	if(vtkOutput > 0) {
		writeVTK(-2);
	}
}

void GaussSeidel2D::initialConditions(int vtkOutput) {
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

	if(vtkOutput > 0) {
		writeVTK(-1);
	}
}

void GaussSeidel2D::printInfo(int iterations, double timeSec) const {
	cout << "runtime: " << timeSec << endl;
	cout << "number of iterations: " << iterations << endl;
	double updatesPerSecond = iterations / timeSec;
	long int gridpoints = ((_nx-1)*(_ny-1));
	double millionGridpointUpdatesPerSecond = gridpoints * updatesPerSecond * 1.0e-6;
	cout << "iterations per second: " << updatesPerSecond << endl;
	cout << "million gridpoint updates per second: " << millionGridpointUpdatesPerSecond << endl;
	long int bytes = (_nx * _ny) * sizeof(double);
	cout << "total data size in  Bytes: " << bytes << endl;
	cout << "total data size in KBytes: " << bytes /1000 << endl;
	cout << "total data size in MBytes: " << bytes /1000000<< endl;
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



