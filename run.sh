#!/bin/bash

nx=20
ny=20
solver=0
numIterations=10
vtkOutputFrequency=5

./another-gauss-seidel $nx $ny $solver $numIterations $vtkOutputFrequency
