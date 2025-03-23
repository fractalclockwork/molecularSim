// File:  cudaSimKernels.h
// -----------------
// This file contains the definition of key cuda kernels required in
// molecular simulation.

#ifndef _cudaSimKernels_h_
#define _cudaSimKernels_h_

const int numThreadsInBlock = 32;

__global__ void cudaReduce(double *sum, double *array);

__global__ void cudaReduceThread(double *sum, double *array);

// MC kernel to calculate LJ energies
__global__
void energyLJKernel(int num, int index, int *type, double length, double *rCut,
     double *sigma, double *epsilon, double *ri, double *r, double *potEnergy);
__global__
void setLJEnergyKernel(int num, int *type, double length, double *rCut,
     double *sigma, double *epsilon, double *r, double *potEnergy);

#endif

