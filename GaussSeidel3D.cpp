/*
 * GaussSeidel3D.cpp
 *
 *  Created on: Feb 23, 2020
 *      Author: tchipevn
 */

#include "GaussSeidel3D.h"
#include <iostream>
#include <fstream>   /* file io */
#include <stdlib.h>  /* srand */
#include <chrono>    /* steady_clock */

using namespace std;
using namespace std::chrono;

void GaussSeidel3D::run(int whichSolver, int numIterations, int vtkOutput) {
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

void GaussSeidel3D::boundaryConditions(int vtkOutput) {
	// left boundary 0, right boundary 1
	for (int z = 0; z < _nz; ++z) {
		for (int y = 0; y < _ny; ++y) {
			val(0,y,z) = 0.0;
			val(_nx-1,y,z) = 1.0;
		}
	}

	// bottom, top boundaries: linear interpolation
	double nx_inverse = 1.0 / (_nx-1);
	for (int z = 0; z < _nz; z+= _nz-1) {

		// only two iterations in z loop

		for (int y = 0; y < _ny; ++y) {
			for (int x = 1; x < _nx-1; ++x) {
				val(x,y,z) = x * nx_inverse;
			}
		}
	}

	// front, back boundaries: linear interpolation
	for (int z = 0; z < _nz; ++z) {
		for (int y = 0; y < _ny; y+=_ny-1) {
			// only two iterations in z loop
			for (int x = 1; x < _nx-1; ++x) {
				val(x,y,z) = x * nx_inverse;
			}
		}
	}

	if(vtkOutput > 0) {
		writeVTK(-2);
	}
}

void GaussSeidel3D::initialConditions(int vtkOutput) {
}

void GaussSeidel3D::printInfo(int iterations, double timeSec) const {
}

void GaussSeidel3D::writeVTK(int iteration) {
	ofstream myfile;
	string filename = "vtk/val_"+ to_string(iteration) + ".vtk";
	myfile.open (filename);
	myfile << "# vtk DataFile Version 2.0\n";
	myfile << "volume example " << endl;
	myfile << "ASCII " << endl;
	myfile << "DATASET STRUCTURED_POINTS " << endl;
	myfile << "DIMENSIONS " << _nx << " " << _ny << " " << _nz << "\n";
	myfile << "ASPECT_RATIO 1 1 1 " << endl;
	myfile << "ORIGIN 0 0 0 " << endl;
	myfile << "POINT_DATA " << _nx * _ny * _nz << "\n";
	myfile << "SCALARS volume_scalars double " << endl;
	myfile << "LOOKUP_TABLE default " << endl;

	for (int z = 0; z < _nz; ++z) {
		for (int y = 0; y < _ny; ++y) {
			for (int x = 0; x < _nx; ++x) {
				myfile << val(x,y,z) << " ";
			}
			myfile << endl;
		}
	}

	myfile.close();
}
