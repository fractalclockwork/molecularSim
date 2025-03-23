// File: nveMC.cpp
// ---------------- 
// File containing functions to implement the
// NVEensemble class.
 
// This code was specifically developed to illustrate concepts in the accompanying book:
// R. J. Sadus, "Molecular Simulation of Fluids: Theory, Algorithms, Object-Orientation,
// and Parallel Computing," 2nd Ed. (Elsevier, Amsterdam, 2023). It can be used freely for
// any not for profit purpose or academic research application. The code has been validated,
// but it would be nonetheless prudent to test it further before publishing any results. 
// Check the book's website for any subsequent updates.

#include <math.h>
#include <string.h>
#include "nveMC.h"
#include "auxfunc.h"
#include "ljenergyMC.h"

// Method: getEnergyLRC
// Usage: n = getEnergyLRC();
// -------------------------- 
// Gets the long range energy correction.

double NVEensembleMC::getEnergyLRC()
{
  return energyLRC;
}

// Method: initialCoord
// Usage: initialCoord();
// ---------------------- 
// Places atoms upon a lattice to get initial coordinates

void NVEensembleMC::initialCoord()
{
  int i, j, x, y, z, offset;
  double cells, dcells, cellL, halfCellL, *tempPos, *position;
 
  // Determine the number of unit cells in each coordinate
  // direction
  dcells = pow(0.25 * (double)numAtom, 1.0/3.0);
  cells =  (int) nearestInt(dcells, 1.0);
 
  // check if numAtom is an non-fcc number of molecules
  // and increase the number of cells if necessary
 
  while((4 * cells * cells * cells) < numAtom)
                 cells = cells + 1;
 
  // Determine length of the unit cell
  cellL = boxLen/ (double) cells;
  halfCellL = 0.5 * cellL;

  // Construct the unit cell
  // point to atoms position
  position = atoms[0]->getPosition();
  position[0] = 0.0;
  position[1] = 0.0;
  position[2]= 0.0;

  position = atoms[1]->getPosition();
  position[0] = halfCellL;
  position[1] = halfCellL;
  position[2]= 0.0;

  position = atoms[2]->getPosition();
  position[0] = 0.0;
  position[1] = halfCellL;
  position[2]= halfCellL;

  position = atoms[3]->getPosition();
  position[0] = halfCellL;
  position[1] = 0.0;
  position[2]= halfCellL;
 
  for(i = 4; i < numAtom; i++)
  {
    position = atoms[i]->getPosition();
    position[0] = 0.0;
    position[1] = 0.0;
    position[2] = 0.0;
  } //init all other atoms to 0
 
  // Build the lattice from the unit cell by
  // repeatly translating the four vectors of
  // the unit cell through a distance cellL in
  // the x, y and z directions
 
  offset = 0;
 
  for (z = 1; z <= cells; z++)
    for (y = 1; y <= cells; y++)
      for (x = 1; x <= cells; x++){
        for (i = 0; i < 4; i++){
          j = i + offset;
          if(j < numAtom){
            tempPos = atoms[j]->getPosition();
            position = atoms[i]->getPosition();
            tempPos[0] = position[0] + cellL * (x-1);
            tempPos[1] = position[1] + cellL * (y-1);
            tempPos[2] = position[2] + cellL * (z-1);
          }
        }
 
      offset = offset + 4;
  }

  // Shift centre of box to the origin.
  for (i = 0; i < numAtom; i++){
      tempPos = atoms[i]->getPosition();
	 tempPos[0] -= halfCellL;
	 tempPos[1] -= halfCellL;
	 tempPos[2] -= halfCellL;
  }
}

// Method: reSetEnergy
// Usage:  reSetEnergy(energy);
// ----------------------------
// Re-sets the energy of the ensemble following
// an acceted move
  
void NVEensembleMC::reSetEnergy(double eng)
{
 potEnergy = eng;
}

// Method: setEnergy
// Usage: setEnergy();
// -------------------
// Set the energy calculations using LJenergy.
// This is an N*N calculation to set the energy
// at the start of the simulation.

void NVEensembleMC::setEnergy()
{
  theEnergy->setEnergy(numAtom, boxLen,
             &potEnergy);
}

// Method: getTrialPotE
// Usage: n = getTrialPotE(index, atoms)
// ---------------------------------
// Performs N calculations to update energy before and after
// index atom is moved.
  
double NVEensembleMC::getTrialPotE(int index, AtomMC **a)
{
  return theEnergy->getTrialPotE(numAtom, boxLen, index, a);
}  

// Method: lrc
// Usage: lrc();
//  ------------
// initiate the force calculations for long range corrections

void NVEensembleMC::lrc()
{
  theEnergy->lrc(numComp, comp, boxVol, &energyLRC);
}

// Method: NVEensemble
// Usage: NVEensemble;
// ------------------- 
//  instantiate the NVEensemble class

NVEensembleMC::NVEensembleMC(AtomMC **atoms, int numComp, int numAtom,
   double totEfixed, double density, 
   double *molFract, int *comp)
   : EnsembleMC(atoms, numComp, numAtom, totEfixed,
     density, molFract, comp)
{
  setVolume();
  setLength();
  setComp();
  initialCoord();
  theEnergy = new LJenergy(atoms);
}
