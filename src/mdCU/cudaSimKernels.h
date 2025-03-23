// File:  cudaSimKernels.h
// -----------------
// This file contains the definition of key cuda kernels required in
// molecular simulation.

#ifndef _cudaSimKernels_h_
#define _cudaSimKernels_h_

// Specify the number of threads to be launched
const int numThreadsInBlock = 32;

// MC kernels to calculate LJ energies
__global__ void energyLJKernel(int num, int index, int* type, double length, double* rCut,
                               double* sigma, double* epsilon, double* ri, double* r,
                               double* potEnergy);
__global__ void setLJEnergyKernel(int num, int* type, double length, double* rCut, double* sigma,
                                  double* epsilon, double* r, double* potEnergy);

// MD kernel to calculate LJ forces
__global__ void setLJForceKernel(int num, int* type, double length, double* rCut, double* sigma,
                                 double* epsilon, double* r, double* force, double* potEnergy,
                                 double* virial);

#endif
