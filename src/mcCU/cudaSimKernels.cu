// File:  cudaSimKernals.cu
// -----------------
// This file contains the implementation details for the CUDA
// kernals required for molecular simulation using GPUs.

// This code was specifically developed to illustrate concepts in the accompanying book:
// R. J. Sadus, "Molecular Simulation of Fluids: Theory, Algorithms, Object-Orientation,
// and Parallel Computing," 2nd Ed. (Elsevier, Amsterdam, 2023). It can be used freely for
// any not for profit purpose or academic research application. The code has been validated,
// but it would be nonetheless prudent to test it further before publishing any results.
// Check the book's website for any subsequent updates.

#include "cudaSimKernels.h"

// Kernel to determine the LJ potential energy experienced by an index atom due
// to interations with all other distinct pairs of atoms. This kernel is used in
// MC simulations. It can be used for either mixtures or pure fluids.

__global__ void energyLJKernel(int num, int index, int* type, double length, double* rCut,
                               double* sigma, double* epsilon, double* ri, double* r,
                               double* potEnergy) {
    int k, nInt, kindi, kindj;
    double rXij, rYij, rZij;        // pair seperation vectors
    double rijSq, rCutSq;           // interatomic seperation squared
    double sigmaSq;                 // squared sigma values
    double sigma2, sigma6, sigma12; // multiples of LJ sigma
    double pot, potE;               // potential between 2 atoms
    double div, cut, eps, sig;

    __shared__ double kernelCache[numThreadsInBlock];

    int tID = threadIdx.x + blockIdx.x * blockDim.x;
    int cacheIndex = threadIdx.x;
    int step = blockDim.x * gridDim.x;

    // Calculate energy experienced by the sellected atom (index)
    // with all other distinct pairs

    kindi = type[index];

    // Loop over all the other atoms
    potE = 0.0;
    for (int j = tID; j < num; j += step) {
        k = 3 * j;

        // Exclude self-interaction
        if (j != index) {
            rXij = ri[0] - r[k];
            rYij = ri[1] - r[k + 1];
            rZij = ri[2] - r[k + 2];

            // Apply periodic boundary conditions
            div = rXij / length;
            if (rXij >= 0)
                nInt = (int)(div + 0.5);
            else
                nInt = -(int)(0.5 - div);

            rXij -= length * nInt;

            div = rYij / length;
            if (rYij >= 0)
                nInt = (int)(div + 0.5);
            else
                nInt = -(int)(0.5 - div);

            rYij -= length * nInt;

            div = rZij / length;
            if (rZij >= 0)
                nInt = (int)(div + 0.5);
            else
                nInt = -(int)(0.5 - div);

            rZij -= length * nInt;

            rijSq = rXij * rXij + rYij * rYij + rZij * rZij;

            kindj = type[j];

            // Account for dissimilar atoms
            if (kindi == kindj)
                rCutSq = rCut[index] * rCut[index];
            else {
                cut = (rCut[index] + rCut[j]) / 2;
                rCutSq = cut * cut;
            }

            // Calculate forces for seprations
            // below the cutoff
            if (rijSq <= rCutSq) {
                // Accound for dissimilar atoms
                if (kindi == kindj) {
                    eps = epsilon[index];
                    sig = sigma[index];
                } else {
                    eps = sqrt(epsilon[index] * epsilon[j]);
                    sig = (sigma[index] + sigma[j]) / 2;
                }

                // Determine potential between i and j
                sigmaSq = sig * sig;
                sigma2 = sigmaSq / rijSq;
                sigma6 = sigma2 * sigma2 * sigma2;
                sigma12 = sigma6 * sigma6;
                pot = sigma12 - sigma6;
                potE += 4 * eps * pot;
            }
        }
    }
    // Set the kernel cache values
    kernelCache[cacheIndex] = potE;

    // Synchronize thread in this block
    __syncthreads();

    // Perform reductions, assuming numThreadsInBox is a power of 2.

    int i = blockDim.x / 2;
    while (i != 0) {
        if (cacheIndex < i)
            kernelCache[cacheIndex] += kernelCache[cacheIndex + i];

        __syncthreads();

        i /= 2;
    }

    if (cacheIndex == 0)
        potEnergy[blockIdx.x] = kernelCache[0];
}

// Kernel to determine the complete LJ energy for all pairs of atoms. This kernel
// is used for MC simulations using CPUs. It can be used for either mixtures or
// pure fluids..

__global__ void setLJEnergyKernel(int num, int* type, double length, double* rCut, double* sigma,
                                  double* epsilon, double* r, double* energy) {
    int i, j, kindi, kindj, nInt;
    double rXi, rYi, rZi;           // position vectors for atom i
    double rXj, rYj, rZj;           // position vectors for atom j
    double rXij, rYij, rZij;        // pair seperation vectors
    double rijSq, rCutSq;           // interatomic seperation squared
    double sigmaSq;                 // squared sigma value
    double sigma2, sigma6, sigma12; // multiples of LJ sigma
    double pot, potE;               // potential between 2 atoms
    double div, cut, sig, eps;

    __shared__ double kernelECache[numThreadsInBlock];

    int tID = threadIdx.x + blockIdx.x * blockDim.x;
    int cacheIndex = threadIdx.x;
    int step = blockDim.x * gridDim.x;

    // perform full N^2 calculation at the begining of the
    // simulation.

    // Calculate energy experienced by atoms due
    // interatomic pair interactions.

    potE = 0.0;

    for (i = tID; i < num; i += step) {
        // Component type
        kindi = type[i];

        // Select position vectors of atom i
        rXi = r[3 * i];
        rYi = r[3 * i + 1];
        rZi = r[3 * i + 2];

        for (j = 0; j < num; j++) {
            if (j != i) {
                // Component type
                kindj = type[j];

                // Select position vectors
                rXj = r[3 * j];
                rYj = r[3 * j + 1];
                rZj = r[3 * j + 2];

                // Calculate pair seperation
                rXij = rXi - rXj;
                rYij = rYi - rYj;
                rZij = rZi - rZj;

                // Apply periodic boundary conditions
                div = rXij / length;
                if (rXij >= 0)
                    nInt = (int)(div + 0.5);
                else
                    nInt = -(int)(0.5 - div);

                rXij -= length * nInt;

                div = rYij / length;
                if (rYij >= 0)
                    nInt = (int)(div + 0.5);
                else
                    nInt = -(int)(0.5 - div);

                rYij -= length * nInt;

                div = rZij / length;
                if (rZij >= 0)
                    nInt = (int)(div + 0.5);
                else
                    nInt = -(int)(0.5 - div);

                rZij -= length * nInt;

                rijSq = rXij * rXij + rYij * rYij + rZij * rZij;

                // Account for dissimilar atoms

                if (kindi == kindj)
                    rCutSq = rCut[i] * rCut[i];
                else {
                    cut = (rCut[i] + rCut[j]) / 2;
                    rCutSq = cut * cut;
                }

                // Calculate forces for seprations
                // below the cutoff

                if (rijSq <= rCutSq) {

                    // Account for dissimilar atoms
                    if (kindi == kindj) {
                        eps = epsilon[i];
                        sig = sigma[i];
                    } else {
                        eps = sqrt(epsilon[i] * epsilon[j]);
                        sig = (sigma[i] + sigma[j]) / 2;
                    }

                    // Determine potential between i and j
                    sigmaSq = sig * sig;
                    sigma2 = sigmaSq / rijSq;
                    sigma6 = sigma2 * sigma2 * sigma2;
                    sigma12 = sigma6 * sigma6;
                    pot = sigma12 - sigma6;
                    potE += 4 * eps * pot;
                }
            }
        }
    }

    // Set the kernel cache values
    kernelECache[cacheIndex] = potE / 2;

    // Synchronize thread in this block
    __syncthreads();

    // Perform reductions, assuming numThreadsInBox is a power of 2.

    int k = blockDim.x / 2;
    while (k != 0) {
        if (cacheIndex < k)
            kernelECache[cacheIndex] += kernelECache[cacheIndex + k];

        __syncthreads();

        k /= 2;
    }

    if (cacheIndex == 0)
        energy[blockIdx.x] = kernelECache[0];
}

__global__ void cudaReduce(double* sum, double* array) {
    extern __shared__ double temp[];

    int tID = threadIdx.x;

    temp[tID] = array[tID + blockIdx.x * blockDim.x];

    for (int i = blockDim.x >> 1; i >= 1; i >>= 1) {
        __syncthreads();

        if (tID < i)
            temp[tID] += temp[tID + i];
    }

    if (tID == 0)
        sum[blockIdx.x] = temp[0];
}

// Kernal to reduce threads is adapted from M. Harris (NVIDIA)
__global__ void cudaReduceThread(double* globalOut, double* globalIn) {
    extern __shared__ double sharedData[];

    // for each thread, load one element from global to shared memory
    unsigned int tID = threadIdx.x;
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    sharedData[tID] = globalIn[i];
    __syncthreads();

    // perform reduction on shared memory assuming threads are powers of 2
    for (unsigned int j = 1; j < blockDim.x; j *= 2) {
        if (tID % (2 * j) == 0)
            sharedData[tID] += sharedData[tID + j];

        __syncthreads();
    }

    // write result for the block to global memory
    if (tID == 0)
        globalOut[blockIdx.x] = sharedData[0];
}
