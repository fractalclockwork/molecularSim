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

#include "ljenergyMC.h"
#include "auxfunc.h"

// Method: setEnergy
// Usage: setEnergy();
// -------------------- 
// The ljenergy method calculates the energy experieced
// by all num atoms due to pairwise interatomic interaction.
// The forces are calculated using the Lennard-Jones (6-12)
// potential. The potential is truncated at a distance
// rCut and long range corrections must be applied outside
// the method to obtain the full contributions to the 
// energy.
 
void LJenergy::setEnergy(int num, double length, 
                double *potEnergy)
{
  int i, j, kindi, kindj;
  double rXi, rYi, rZi;           // position vectors for atom i
  double rXj, rYj, rZj;           // position vectors for atom j
  double rXij, rYij, rZij;        // pair seperation vectors
  double rijSq;                   // interatomic seperation squared
  double rCutSq, sigmaSq;         // squared cut and sigma values
  double sigma2, sigma6, sigma12; // multiples of LJ sigma
  double pot , eps;               // potential between 2 atoms
  double *tempPos, *tempPosj;

  *potEnergy = 0.0;

  // perform full N^2 calculation at the begining of the
  // simulation.
  
  // Calculate energy experienced by atoms due
  // interatomic pair interactions.
  for(i = 0; i < num - 1; i++)
  {
    // Select molecule i
    kindi = atoms[i]->getType();
    tempPos = atoms[i]->getPosition();

    rXi   = tempPos[0];
    rYi   = tempPos[1];
    rZi   = tempPos[2];

    for(j = i + 1; j < num; j++)
    {
      kindj = atoms[j]->getType();
      tempPosj = atoms[j]->getPosition();

      rXj   = tempPosj[0];
      rYj   = tempPosj[1];
      rZj   = tempPosj[2];

      // Calculate pair seperation
      rXij = rXi - rXj;
      rYij = rYi - rYj;
      rZij = rZi - rZj;

      // Apply periodic boundary
      rXij -= length * nearestInt(rXij, length);
      rYij -= length * nearestInt(rYij, length);
      rZij -= length * nearestInt(rZij, length);

      rijSq = rXij * rXij + rYij * rYij + rZij * rZij;
      rCutSq = rCut[kindi][kindj] * rCut[kindi][kindj];

      // Calculate energies for seprations
      // below the cutoff
      if (rijSq < rCutSq)
      {
  	    // Determine potential between i and j
	    sigmaSq = sigma[kindi][kindj] * sigma[kindi][kindj];
	    sigma2  = sigmaSq/rijSq;
	    sigma6  = sigma2 * sigma2 * sigma2;
	    sigma12 = sigma6 * sigma6;
	    pot     = sigma12 - sigma6;
        eps     = epsilon[kindi][kindj];
	    *potEnergy   += eps * pot;
      }
    }
  }
  
  *potEnergy *= 4;
}

// Method: getTrialPotE
// Usage: getTrialPotE();
// -------------------- 
// The ljenergy method calculates the energy experieced
// by a given atom (index) due to pairwise interatomic interaction.
// The energies are calculated using the Lennard-Jones (6-12)
// potential. The potential is truncated at a distance
// rCut and long range corrections must be applied outside
// the method to obtain the full contributions to the 
// energy.
 
double LJenergy::getTrialPotE(int num, double length, int index, AtomMC **atom)
{
  int i, j, kindi, kindj;
  double rXi, rYi, rZi;           // position vectors for atom i
  double rXj, rYj, rZj;           // position vectors for atom j
  double rXij, rYij, rZij;        // pair seperation vectors
  double rijSq;                   // interatomic seperation squared
  double rCutSq, sigmaSq;         // squared cut and sigma values
  double sigma2, sigma6, sigma12; // multiples of LJ sigma
  double pot , potE, eps;         // potential between 2 atoms
  double *tempPos, *tempPosj;

  potE = 0.0;

  // Calculate energy experienced by the sellected atom (index)
  // with all other distinct pairs
    
  // Properties of index atom.
  kindi = atom[index]->getType();
  tempPos = atom[index]->getTrialPosition();

  rXi   = tempPos[0];
  rYi   = tempPos[1];
  rZi   = tempPos[2];

  // Loop over all the other atoms
  for(j = 0; j < num; j++)
  {
    //Exclude self-interaction
    if(j != index)
    {
      kindj = atom[j]->getType();
      tempPosj = atom[j]->getPosition();

      rXj   = tempPosj[0];
      rYj   = tempPosj[1];
      rZj   = tempPosj[2];

      // Calculate pair seperation
      rXij = rXi - rXj;
      rYij = rYi - rYj;
      rZij = rZi - rZj;

      // Apply periodic boundary
      rXij -= length * nearestInt(rXij, length);
      rYij -= length * nearestInt(rYij, length);
      rZij -= length * nearestInt(rZij, length);

      rijSq = rXij * rXij + rYij * rYij + rZij * rZij;
      rCutSq = rCut[kindi][kindj] * rCut[kindi][kindj];

      // Calculate energies for seprations
      // below the cutoff
      if (rijSq < rCutSq)
      {
  	   // Determine potential between i and j
	   sigmaSq = sigma[kindi][kindj] * sigma[kindi][kindj];
	   sigma2  = sigmaSq/rijSq;
	   sigma6  = sigma2 * sigma2 * sigma2;
	   sigma12 = sigma6 * sigma6;
	   pot     = sigma12 - sigma6;
       eps     = epsilon[kindi][kindj];
	   potE    += eps * pot;
      }
    }
  }

  return 4*potE;
    
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
LJenergy::LJenergy(AtomMC **theAtoms)
:Energy(theAtoms)
{
  epsilon = atoms[0]->getEpsilon();
  sigma   = atoms[0]->getSigma();
  rCut    = atoms[0]->getrCutOff();
}
