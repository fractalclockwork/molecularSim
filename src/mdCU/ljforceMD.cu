// File: LJforceMD.cu
// -------------------
// File containing methods to implement the LJforce
// class.

// This code was specifically developed to illustrate concepts in the accompanying book:
// R. J. Sadus, "Molecular Simulation of Fluids: Theory, Algorithms, Object-Orientation,
// and Parallel Computing," 2nd Ed. (Elsevier, Amsterdam, 2023). It can be used freely for
// any not for profit purpose or academic research application. The code has been validated,
// but it would be nonetheless prudent to test it further before publishing any results.
// Check the book's website for any subsequent updates.

#include "cudaSimKernels.h"
#include "ljforceMD.h"
#include <iostream>
#include <stdlib.h>

using namespace std;

// Method: setForce
// Usage: setForce();
// ------------------
// The ljForce method uses a GPU kernel to calculate the force
// experienced by all num atoms due to pairwise interatomic interaction.
// The forces are calculated using the kernel for the Lennard-Jones (6-12)
// potential. This can be replaced by kernels for other potentials.
// The potential is truncated at a distance of rCut and
// long range corrections must be applied outside of the
// the method to obtain the full contributions to both
// energy and the virial.

void LJforce::setForce(int num, double length, double* potEnergy, double* virial) {

    int* kind_dev;               // atom type;
    double* r_dev;               // atom position vectors
    double *energy, *energy_dev; // potential energy;
    double *vir, *vir_dev;       // virial
    double* force_dev;           // forces
    double *sigma_dev, *epsilon_dev, *rCut_dev;
    double potE, virE;

    // Distribute the calculations among the threads
    int blockSize = numThreadsInBlock;
    int numBlocks = (num + blockSize - 1) / blockSize;

    // Reinitalise force array prior to kernel execution
    for (int i = 0; i < num; i++) {
        atoms[0]->f[3 * i] = 0;
        atoms[0]->f[3 * i + 1] = 0;
        atoms[0]->f[3 * i + 2] = 0;
    }

    // Allocate memory on the CPU side
    energy = (double*)malloc(numBlocks * sizeof(double));
    vir = (double*)malloc(numBlocks * sizeof(double));

    // Allocate memory on the GPU
    cudaMalloc((void**)&r_dev, 3 * num * sizeof(double));
    cudaMalloc((void**)&force_dev, 3 * num * sizeof(double));
    cudaMalloc((void**)&energy_dev, numBlocks * sizeof(double));
    cudaMalloc((void**)&vir_dev, numBlocks * sizeof(double));
    cudaMalloc((void**)&kind_dev, num * sizeof(int));
    cudaMalloc((void**)&sigma_dev, num * sizeof(double));
    cudaMalloc((void**)&epsilon_dev, num * sizeof(double));
    cudaMalloc((void**)&rCut_dev, num * sizeof(double));

    // Copy arrays to the GPU
    cudaMemcpy(r_dev, atoms[0]->r, 3 * num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(force_dev, atoms[0]->f, 3 * num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(kind_dev, atoms[0]->kind, num * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(sigma_dev, atoms[0]->sigii, num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(epsilon_dev, atoms[0]->epsii, num * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(rCut_dev, atoms[0]->rCutii, num * sizeof(double), cudaMemcpyHostToDevice);

    // Execute the kernel on the GPU
    setLJForceKernel<<<numBlocks, blockSize>>>(num, kind_dev, length, rCut_dev, sigma_dev,
                                               epsilon_dev, r_dev, force_dev, energy_dev, vir_dev);

    // Copy the 'energy' 'virial' and 'force' arrays back from the GPU to the CPU
    cudaMemcpy(energy, energy_dev, numBlocks * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(vir, vir_dev, numBlocks * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(atoms[0]->f, force_dev, 3 * num * sizeof(double), cudaMemcpyDeviceToHost);

    // Sum energy contributions of the CPU side
    potE = 0;
    virE = 0;
    for (int k = 0; k < numBlocks; k++) {
        potE += energy[k];
        virE += vir[k];
    }

    *potEnergy = potE;
    *virial = virE;

    // Free memory on the GPU side
    cudaFree(energy_dev);
    cudaFree(r_dev);
    cudaFree(force_dev);
    cudaFree(kind_dev);
    cudaFree(sigma_dev);
    cudaFree(epsilon_dev);
    cudaFree(rCut_dev);

    // Free memory of the CPU side,
    free(energy);
    free(vir);
}

// Method: lrc
// Usage: n = lrc(num, nComp, volume, energyLRC, virialLRC);
// ---------------------------------------------------------
// ljLRC calculates the long range correction terms to both
// the energy and virial for an ensemble interacting via the
// lennard-Jones (12-6) potential.

void LJforce::lrc(int num, int* nComp, double boxV, double* energyLRC, double* virialLRC) {
    int i, j;
    const double PI = 3.141592654;
    double sigmaSq, sig3, rCutSq, rCut3, sig3R, sig9R, den;

    *energyLRC = 0.0;
    *virialLRC = 0.0;

    for (i = 0; i < num; i++)
        for (j = 0; j < num; j++) {
            sigmaSq = sigma[i][j] * sigma[i][j];
            rCutSq = rCut[i][j] * rCut[i][j];
            sig3 = sigma[i][j] * sigmaSq;
            rCut3 = rCut[i][j] * rCutSq;
            sig3R = sig3 / rCut3;
            sig9R = sig3R * sig3R * sig3R;
            den = nComp[i] * nComp[j] * sig3 / boxV;
            *energyLRC += (8.0 / 9.0) * den * PI * epsilon[i][j] * (sig9R - 3 * sig3R);
            *virialLRC += (16.0 / 9.0) * den * PI * epsilon[i][j] * (2 * sig9R - 3 * sig3R);
        }
}

// constructor
LJforce::LJforce(Atom** theAtoms) : Force(theAtoms) {
    epsilon = atoms[0]->getEpsilon();
    sigma = atoms[0]->getSigma();
    rCut = atoms[0]->getrCutOff();
}
