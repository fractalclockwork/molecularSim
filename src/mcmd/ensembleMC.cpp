// File: ensembleMC.cpp
// -------------------- 
// File containing functions to implement the
// Abstract Ensemble class.

// This code was specifically developed to illustrate concepts in the accompanying book:
// R. J. Sadus, "Molecular Simulation of Fluids: Theory, Algorithms, Object-Orientation,
// and Parallel Computing," 2nd Ed. (Elsevier, Amsterdam, 2023). It can be used freely for
// any not for profit purpose or academic research application. The code has been validated,
// but it would be nonetheless prudent to test it further before publishing any results. 
// Check the book's website for any subsequent updates.

#include "ensembleMC.h"
#include <math.h>

// setMethods
// ----------
// These methods assign values to the volume, length, and component number
// for the ensemble.  They do not take any parameters, and are called
// during the construction phase of the ensemble

// Method: setVolume
// Usage: setVolume();
// ------------------- 
//  Determines the volume of the simulation box.

void EnsembleMC::setVolume()
{
  boxVol = numAtom/density;
}

// Method: setLength
// Usuage: setLength();
// -------------------- 
// Determines the length of the simulation box.

void EnsembleMC::setLength()
{
  boxLen = pow(boxVol, 1.0/3.0);
}

// Method: setComp();
// Usage: setComp();
// ------------------
//  Determines the number of atoms of each type.

void EnsembleMC::setComp()
{
  int i, sum = 0;

  for(i =0; i < numComp - 1; i++){
     comp[i] = (int)(molFract[i] * numAtom);
		sum += comp[i];
  }
  comp[numComp - 1] = numAtom - sum;
}

void EnsembleMC::setKineticE(double kNew)
{
 kineticE = kNew;
}

// getMethods
// -----------
// These methods return the numebr of atoms, an array describing the
// number of each atoms of each type, the number of components, the
// temperature of the ensemble, the volume of the ensemble, the length of
// the ensemble, the potential energy, the virial, and the kinetic energy
// of the ensemble.

// Method: getNumAtom
// Usage: n = getNumAtom();
// ------------------------ 
// Gets the total number of molecules in the ensemble.

int EnsembleMC::getNumAtom()
{
  return numAtom;
}

// Function: getNumComp
// Usage: n = getNumComp();
// ------------------------_ 
// Gets the number of components in the ensemble.

int EnsembleMC::getNumComp()
{
  return numComp;
}

// Method: getComp
// Usage: p = getComp();
// ---------------------

// Gets the number of each individual component.

int * EnsembleMC::getComp()
{
  return comp;
}

// Method: getVolume
// Usage: getVolume();
// -------------------
// Gets the volume of the simulation box.

double EnsembleMC::getVolume()
{
  return boxVol;
}

// Method: getLength
// Usage: getLength();
// ------------------- 
// Gets the length of the simulation box.

double EnsembleMC::getLength()
{
  return boxLen;
}

// Method: getPotEnergy
// Usage: getPotEnergy();
// ---------------------- 
// Gets potential energy following force calculation.

double EnsembleMC::getPotEnergy()
{
  return potEnergy;
}
// Method getKineticE
// Usage n = getKineticE();
// ------------------------ 
// Gets the kinetic energy.

double EnsembleMC::getKineticE()
{
  return kineticE;
}

// Method: gettotEfixed
// Usage: n = gettotEfixed();
// --------------------------
// Get the pre-determined total energy of the
// ensemble.

double EnsembleMC::gettotEfixed()
{
  return totEfixed;
}

// Method: getAtoms
// Usage: n = getAtoms();
// ---------------------- 
//  return array of atoms

AtomMC **EnsembleMC::getAtoms()
{
  return atoms;
}

// Method: Ensemble
// Usage: Ensemble;
// ---------------- 
// Constructor for the Ensemble class

EnsembleMC::EnsembleMC(AtomMC **theAtoms, int nComp, int nAtoms, 
	double totE, double dens, double *mol, int *cmp)
{
  atoms = theAtoms;
  numComp = nComp;
  numAtom = nAtoms;
  totEfixed = totE;
  density = dens;
  molFract = mol;
  comp = cmp;
  potEnergy = 0.0;
  trialPotE = 0.0;
}

// Destructor for Ensemble class
EnsembleMC::~EnsembleMC()
{
}


