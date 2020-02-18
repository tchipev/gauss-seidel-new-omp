/*
 * GaussSeidel2D.cpp
 *
 *  Created on: Feb 9, 2020
 *      Author: tchipevn
 */

#include "GaussSeidel2D.h"
#include <iostream>
#include <fstream>

using namespace std;

void GaussSeidel2D::run(int whichSolver, int numIterations, int vtkOutput) {
	std::cout << "hi from " << whichSolver << std::endl;

	for(int y=1; y< _ny-1; ++y) {
		for(int x=1; x < _nx-1; ++x) {
			process9(x,y);
		}
	}

	writeVTK(10);
	writeVTK(11);
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



