// File: LJforceMD.cpp
// ------------------- 
// File containing methods to implement the LJforce
// class.
 
// This code was specifically developed to illustrate concepts in the accompanying book:
// R. J. Sadus, "Molecular Simulation of Fluids: Theory, Algorithms, Object-Orientation,
// and Parallel Computing," 2nd Ed. (Elsevier, Amsterdam, 2023). It can be used freely for
// any not for profit purpose or academic research application. The code has been validated,
// but it would be nonetheless prudent to test it further before publishing any results. 
// Check the book's website for any subsequent updates.
 
#include "ljforceMD.h"
#include "auxfunc.h"
#include <iostream>
#include <fstream>

using namespace std;
 
// Method: setForce
// Usage: setForce();
// ------------------ 
// The ljForce function calculates the force experieced
// by all num atoms due to pairwise interatomic interaction.
// The forces are calculated using the Lennard-Jones (6-12)
// potential. The potential is truncated at a distance
// rCut and long range corrections must be pplied outside
// the function to obtain the full contributions to both
// energy and the virial.
 
void LJforce::setForce(int num, double length,
              double *potEnergy, double *virial)
{
  int i, j, kindi, kindj;
  double rXi, rYi, rZi;           // position vectors for atom i
  double rXj, rYj, rZj;           // position vectors for atom j
  double rXij, rYij, rZij;        // pair seperation vectors
  double rijSq;                   // interatomic seperation squared
  double fXi, fYi, fZi;           // force vectors for atom i
  double fij, fXij, fYij, fZij;   // force between atoms i and j
  double rCutSq, sigmaSq;         // squared cut and sigma values
  double sigma2, sigma6, sigma12; // multiples of LJ sigma
  double pot;                     // Potential between 2 atoms
  double *tempPos;                // postion vectors
  double *tempPosj;               // postion vector
  double potECal, virialCal; 
  double *fX = new double [num]; 
  double *fY = new double [num];
  double *fZ = new double [num];

  // Zero the force  arrays.
  for(i = 0; i < num; i++)
   {
     fX[i] = 0.0;
     fY[i] = 0.0;
     fZ[i] = 0.0;
   }

  potECal = 0.0;
  virialCal = 0.0;

    
  //Calculate forces experienced by atoms due
  //interatomic pair interactions.

  for(i = 0; i < num; i++)
  {
    //Select molecue i
    kindi = atoms[i]->getType();
    tempPos = atoms[i]->getPosition();
    
    rXi   = tempPos[0];
    rYi   = tempPos[1];
    rZi   = tempPos[2];

    fXi = fX[i];
    fYi = fY[i];
    fZi = fZ[i];

    for(j = i + 1; j < num; j++)
    {
       kindj = atoms[j]->getType();
       tempPosj = atoms[j]->getPosition();

       rXj   = tempPosj[0];
       rYj   = tempPosj[1];
       rZj   = tempPosj[2];

       //Calculate pair seperation
       rXij = rXi - rXj;
       rYij = rYi - rYj;
       rZij = rZi - rZj;

       //Apply periodic boundary
       rXij -= length * nearestInt(rXij, length);
       rYij -= length * nearestInt(rYij, length);
       rZij -= length * nearestInt(rZij, length);
        
       rijSq = rXij * rXij + rYij * rYij + rZij * rZij;
       rCutSq = rCut[kindi][kindj] * rCut[kindi][kindj];

       //Calculate forces for seprations
       //below the cutoff
       if (rijSq <= rCutSq)
       {
          //Determine potential between i and j
          sigmaSq = sigma[kindi][kindj] * sigma[kindi][kindj];
          sigma2  = sigmaSq/rijSq;
   	      sigma6  = sigma2 * sigma2 * sigma2;
	      sigma12 = sigma6 * sigma6;
	      pot     = sigma12 - sigma6;

          potECal   += 4*epsilon[kindi][kindj] *pot;
          virialCal += 24*epsilon[kindi][kindj] *(pot + sigma12)/3;

	      //Calculate forces between atoms
          fij  = 24*epsilon[kindi][kindj] * (pot + sigma12)/rijSq;
	      fXij = fij * rXij;
	      fYij = fij * rYij;
	      fZij = fij * rZij;
          fXi += fXij;
	      fYi += fYij;
	      fZi += fZij;
			
          //Apply Newton's third law
          fX[j] -= fXij;
          fY[j] -= fYij;
          fZ[j] -= fZij;
 
       }

     } 

    fX[i]  = fXi;
    fY[i]  = fYi;
    fZ[i]  = fZi;

  } 

  // Store forces
  for(i = 0; i < num; i++)
    atoms[i]->setForce(fX[i], fY[i], fZ[i]);

  *potEnergy = potECal;
  *virial    = virialCal;

  delete [] fX;
  delete [] fY;
  delete [] fZ;

}

// Method: lrc
// Usage: n = lrc(num, nComp, volume, energyLRC, virialLRC);
// --------------------------------------------------------- 
// ljLRC calculates the long range correction terms to both
// the energy and virial for an ensemble interacting via the
// lennard-Jones (12-6) potential.

void LJforce::lrc(int num, int *nComp, double boxV,
		double *energyLRC, double *virialLRC)
{
  int i, j;
  const double PI = 3.141592654;
  double sigmaSq, sig3, rCutSq, rCut3, sig3R, sig9R, den;

  *energyLRC = 0.0;
  *virialLRC = 0.0;

  for(i = 0; i < num; i++)
    for(j = 0; j < num; j++)
    {
      sigmaSq   = sigma[i][j] * sigma[i][j];
      rCutSq    = rCut[i][j] * rCut[i][j];
      sig3      = sigma[i][j] * sigmaSq;
      rCut3     = rCut[i][j] * rCutSq;
      sig3R     = sig3/rCut3;
      sig9R     = sig3R * sig3R * sig3R;
      den       = nComp[i] * nComp[j] * sig3/boxV;
      *energyLRC += (8/9.0) * den * PI * epsilon[i][j]
			* (sig9R - 3 * sig3R);
      *virialLRC += (16/9.0) * den * PI * epsilon[i][j]
			* (2 * sig9R - 3 * sig3R);
   }
}

// constructor
LJforce::LJforce(Atom **theAtoms)
:Force(theAtoms)
{
  epsilon = atoms[0]->getEpsilon();
  sigma   = atoms[0]->getSigma();
  rCut    = atoms[0]->getrCutOff();
}
