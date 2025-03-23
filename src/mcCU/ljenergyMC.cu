// File: ljenergyMC.cpp 
// -------------------- 
// File containing functions to implement the LJenergy
// class.

// This code was specifically developed to illustrate concepts in the accompanying book:
// R. J. Sadus, "Molecular Simulation of Fluids: Theory, Algorithms, Object-Orientation,
// and Parallel Computing," 2nd Ed. (Elsevier, Amsterdam, 2023). It can be used freely for
// any not for profit purpose or academic research application. The code has been validated,
// but it would be nonetheless prudent to test it further before publishing any results.
// Check the book's website for any subsequent updates.

#include "cudaSimKernels.h"
#include "ljenergyMC.h"
#include <iostream>
#include <stdlib.h>

using namespace std;

// Method: getTrialPotE
// Usage: getTrialPotE();
// -------------------- 
// The ljenergy method uses a GPU kernel to calculates the energy experieced
// by a given atom (index) due to pairwise interatomic interaction.
// The energies are calculated using the Lennard-Jones (6-12)
// GPU kernel. This can be replaced by kernels for any other potential
// The potential is truncated at a distance of rCut and 
// long range corrections must be applied outside of the method
// to obtain the full contributions to the energy. 
 
double LJenergy::getTrialPotE(int num, double length, int index, Atom **atom)
{
  int    *kind_dev;                   // atom type
  double *ri, *ri_dev;                // position vecotors of atom i
  double *r_dev;                      // position vectors for all atoms
  double *energy, *energy_dev;        // potential energy; 
  double *position, potE;
  double *sigma_dev, *epsilon_dev;
  double *rCut_dev;

  // Distribute the calculation among the threads
  int blockSize = numThreadsInBlock;
  int numBlocks = (num + blockSize - 1)/blockSize;

  // Allocate memory on the CPU side
  ri = (double *) malloc(3*sizeof(double)); //ri[0], ri[1] and ri[2] are the x, y, z components
  energy = (double *) malloc(numBlocks*sizeof(double));

  // Properties of the index atom
  position  = atom[index]->getTrialPosition();
  ri[0]     = position[0];
  ri[1]     = position[1];
  ri[2]     = position[2];

  // Allocate memory on the GPU
  cudaMalloc((void **) &ri_dev,      3*sizeof(double));
  cudaMalloc((void **) &r_dev,       3*num*sizeof(double));
  cudaMalloc((void **) &energy_dev,  numBlocks*sizeof(double));
  cudaMalloc((void **) &kind_dev,    num*sizeof(int));
  cudaMalloc((void **) &sigma_dev,   num*sizeof(double));
  cudaMalloc((void **) &epsilon_dev, num*sizeof(double));
  cudaMalloc((void **) &rCut_dev,    num*sizeof(double));

  // Copy arrays to the GPU
  cudaMemcpy(ri_dev,     ri,              3*sizeof(double),     cudaMemcpyHostToDevice);
  cudaMemcpy(r_dev,      atoms[0]->r,     3*num*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(kind_dev,   atoms[0]->kind,  num*sizeof(int),      cudaMemcpyHostToDevice);
  cudaMemcpy(sigma_dev,  atoms[0]->sigii, num*sizeof(double),   cudaMemcpyHostToDevice);
  cudaMemcpy(epsilon_dev,atoms[0]->epsii, num*sizeof(double),   cudaMemcpyHostToDevice);
  cudaMemcpy(rCut_dev,   atoms[0]->rCutii,num*sizeof(double),   cudaMemcpyHostToDevice);
 
  // Execute  kernel on the GPU to determine the energy of the index atom with all
  // other distinct pairs
  energyLJKernel <<<numBlocks, blockSize>>> (num, index, kind_dev, length, rCut_dev, 
                               sigma_dev, epsilon_dev, ri_dev, r_dev, energy_dev); 

  // Copy the 'energy' array back from the GPU to the CPU
  cudaMemcpy(energy, energy_dev, numBlocks*sizeof(double), cudaMemcpyDeviceToHost);

  // Sum energy contributions of the CPU side
  potE = 0;
  for(int k = 0; k < numBlocks; k++)
    potE +=  energy[k]; 

  // Free memory on the GPU side
  cudaFree(energy_dev);
  cudaFree(ri_dev);
  cudaFree(r_dev);
  cudaFree(kind_dev);
  cudaFree(sigma_dev);
  cudaFree(epsilon_dev);
  cudaFree(rCut_dev);

  // Free memory of the CPU side
  free(energy);
  free(ri);

  return potE;
    
}

// Method: setEnergy
// Usage: setEnergy();
// -------------------- 
// The ljenergy method  uses a GPU kernel to calculates the energy experieced
// by all num atoms due to pairwise interatomic interaction.
// The forces are calculated using thekernel for the Lennard-Jones (6-12)
// This can be replaced by kernels for any other potential.
// The potential is truncated at a distance of rCut and long range
// corrections must be applied outsid e the method to obtain the
// the full contributions to the energy. 
 
void LJenergy::setEnergy(int num, int numComp, double length, 
                double *potEnergy)
{
  int    *kind_dev;                            //atom type;
  double *r_dev;                               // atom position vectors
  double *energy, *energy_dev;                 // potential energy; 
  double *sigma_dev, *epsilon_dev, *rCut_dev;
  double  potE;

  // Distribute the calculations among the threads  
  int blockSize = numThreadsInBlock;
  int numBlocks = (num + blockSize - 1)/blockSize;

  // Allocate memory on the CPU side
  energy = (double *) malloc(numBlocks*sizeof(double));

  // Allocate memory on the GPU
  cudaMalloc((void **) &r_dev,       3*num*sizeof(double));
  cudaMalloc((void **) &energy_dev,  numBlocks*sizeof(double));
  cudaMalloc((void **) &kind_dev,    num*sizeof(int));
  cudaMalloc((void **) &sigma_dev,   num*sizeof(double));
  cudaMalloc((void **) &epsilon_dev, num*sizeof(double));
  cudaMalloc((void **) &rCut_dev,    num*sizeof(double));

  // Copy arrays to the GPU
  cudaMemcpy(r_dev,      atoms[0]->r,     3*num*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(kind_dev,   atoms[0]->kind,  num*sizeof(int),      cudaMemcpyHostToDevice);
  cudaMemcpy(sigma_dev,  atoms[0]->sigii, num*sizeof(double),   cudaMemcpyHostToDevice);
  cudaMemcpy(epsilon_dev,atoms[0]->epsii, num*sizeof(double),   cudaMemcpyHostToDevice);
  cudaMemcpy(rCut_dev,   atoms[0]->rCutii,num*sizeof(double),   cudaMemcpyHostToDevice);

  // Execute the kernel on the GPU
  setLJEnergyKernel <<<numBlocks, blockSize>>> (num, kind_dev, length, 
                    rCut_dev, sigma_dev, epsilon_dev, r_dev, energy_dev); 

  // Copy the 'energy' array back from the GPU to the CPU
  cudaMemcpy(energy, energy_dev, numBlocks*sizeof(double), cudaMemcpyDeviceToHost);

  // Sum energy contributions of the CPU side
  potE = 0;
  for(int k = 0; k < numBlocks; k++)
    potE +=  energy[k]; 

  *potEnergy = potE;

  // Free memory on the GPU side
  cudaFree(energy_dev);
  cudaFree(r_dev);
  cudaFree(kind_dev);
  cudaFree(sigma_dev);
  cudaFree(epsilon_dev);
  cudaFree(rCut_dev);

  // Free memory of the CPU side
  free(energy);

}

// Method: lrc
// Usage: n = lrc(num, nComp, volume, energyLRC);
// --------------------------------------------------------- 
// ljLRC calculates the long range correction terms to both
// the energy  for an ensemble interacting via the
// lennard-Jones (12-6) potential.

void LJenergy::lrc(int num, int *nComp, double boxV,
		double *energyLRC)
{
  int i, j;
  const double PI = 3.141592654;
  double eps, sigmaSq, sig3, rCutSq, rCut3, sig3R, sig9R, den;

  *energyLRC = 0.0;

  for(i = 0; i < num; i++)
      for(j = 0; j < num; j++){
          eps       = epsilon[i][j];
          sigmaSq   = sigma[i][j] * sigma[i][j];
          rCutSq    = rCut[i][j] * rCut[i][j];
          sig3      = sigma[i][j] * sigmaSq;
          rCut3     = rCut[i][j] * rCutSq;
          sig3R     = sig3/rCut3;
          sig9R     = sig3R * sig3R * sig3R;
          den       = nComp[i] * nComp[j] * sig3/boxV;
          *energyLRC += (8.0/9.0) * den * PI * eps
                    * (sig9R - 3 * sig3R);
      }
}

// constructor
LJenergy::LJenergy(Atom **theAtoms)
:Energy(theAtoms)
{
  epsilon = atoms[0]->getEpsilon();
  sigma   = atoms[0]->getSigma();
  rCut    = atoms[0]->getrCutOff();
}
