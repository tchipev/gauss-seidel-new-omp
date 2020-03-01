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
#include <cassert>

using namespace std;
using namespace std::chrono;

void GaussSeidel2D::run(int whichSolver, int numIterations, int vtkOutput) {

	printSchemeInfo(whichSolver);

	int tenPercentOfIterations = numIterations / 10;
	if (tenPercentOfIterations == 0) {
		tenPercentOfIterations = 1;
	}

	steady_clock::time_point t1 = steady_clock::now();
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
		case 4:
			sumDiff2 = c04_hcpTraversal2();
			break;
		case 5:
			sumDiff2 = c04_hcpTraversal3();
			break;
		case 6:
			sumDiff2 = sli_along_y();
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

void GaussSeidel2D::printSchemeInfo(int solver) const {

	int ompEnvNumThreads;
	#pragma omp parallel
	{
		#pragma omp master
		ompEnvNumThreads = omp_get_num_threads();
	}

	cout << "SOLVER INFO: " << endl;
	switch(solver) {
	case 0:
		cout << "Using basic traversal:" << endl;
		cout << "  iterating over the dimensions in the normal fashion: for(y) {for(x) {}}." << endl;
		cout << "Unsafe OpenMP parallelization!" << endl;
		cout << "Will use: " << ompEnvNumThreads << " threads."<< endl;
		break;
	case 1:
		cout << "Using basic traversal:" << endl;
		cout << "  but iterating over the dimensions in a SLOW fashion: for(x) {for(y) {}}." << endl;
		cout << "Unsafe OpenMP parallelization!" << endl;
		cout << "Will use: " << ompEnvNumThreads << " threads."<< endl;
		break;
	case 2:
		cout << "Using 2D variant of c08." << endl;
		cout << "Will use: " << ompEnvNumThreads << " threads."<< endl;
		break;
	case 3:
		cout << "Using suboptimal 2D variant of c04 (3 colours):" << endl;
		cout << "  only outer loop is parallelized." << endl;
		cout << "Will use: " << ompEnvNumThreads << " threads."<< endl;
		break;
	case 4:
		cout << "Using 2D variant of c04 (3 colours): " << endl;
		cout << "  both loops parallelized. Variant 1. " << endl;
		cout << "Will use: " << ompEnvNumThreads << " threads."<< endl;
		break;
	case 5:
		cout << "Using 2D variant of c04 (3 colours): " << endl;
		cout << "  both loops parallelized. Variant 2. " << endl;
		cout << "Will use: " << ompEnvNumThreads << " threads."<< endl;
		break;
	case 6:
		cout << "Using 2D variant of sli: " << endl;
		cout << "  slicing ONLY along Y dimension." << endl;
		if (_Tx * _Ty != ompEnvNumThreads) {
			cout << "Thread arguments and env var mismatch!" << endl;
			assert(_Tx * _Ty == ompEnvNumThreads);
		}
		cout << "OMP_NUM_THREADS: " << ompEnvNumThreads << endl;
		cout << "Scheme supports maximally " << _slice1dY.getMaxThreads() << " threads." << endl;
		cout << "Will use: " << _slice1dY.getActualThreads() << " threads."<< endl;
		break;
	case 7:
		assert(_Tx == _slice1dX.getActualThreads());
		assert(_Ty == _slice1dY.getActualThreads());
		break;
	}

	cout << endl;
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



