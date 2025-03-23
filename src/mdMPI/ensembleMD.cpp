// File: EnsembleMD.cpp
// -------------------
// File containing methods to implement the
// Abstract Ensemble class.

// This code was specifically developed to illustrate concepts in the accompanying book:
// R. J. Sadus, "Molecular Simulation of Fluids: Theory, Algorithms, Object-Orientation,
// and Parallel Computing," 2nd Ed. (Elsevier, Amsterdam, 2023). It can be used freely for
// any not for profit purpose or academic research application. The code has been validated,
// but it would be nonetheless prudent to test it further before publishing any results.
// Check the book's website for any subsequent updates.

#include "ensembleMD.h"
#include <iostream>
#include <math.h>

using namespace std;

// setMethods
// ----------
// These methods assign values to the volume, length, and component number
// for the ensemble.  They do not take any parameters, and are caled
// during the construction phase of the ensemble

// Method: setVolume
// Usage: setVolume();
// -------------------
// Determines the volume of the simulation box.

void Ensemble::setVolume() { boxVol = numAtom / density; }

// Method: setLength
// Usuage: setLength();
// --------------------
// Determines the length of the simulation box.

void Ensemble::setLength() { boxLen = pow(boxVol, 1.0 / 3); }

// Method: setComp();
// Usage: setComp();
// ------------------
// Determines the number of atoms of each type.

void Ensemble::setComp() {
    int i, sum = 0;

    for (i = 0; i < numComp - 1; i++) {
        comp[i] = (int)(molFract[i] * numAtom);
        sum += comp[i];
    }

    comp[numComp - 1] = numAtom - sum;
}

// getMethods
// ----------
// These methods return the numebr of atoms, an array describing the
// number of each atoms of each type, the number of components, the
// temperature of the ensemble, the volume of the ensemble, the length of
// the ensemble, the potential energy, the virial, and the kinetic energy
// of the ensemble.

// Method: getNumAtom
// Usage: n = getNumAtom();
// ------------------------
// Gets the total number of molecules in the ensemble.

int Ensemble::getNumAtom() { return numAtom; }

// Method: getTemp
// Usage: getTemp();
// -----------------
// Gets the temperature of the ensemble.

double Ensemble::getTemp() { return temperature; }

// Method: getNumComp
// Usage: n = getNumComp();
// ------------------------
// Gets the number of components in the ensemble.

int Ensemble::getNumComp() { return numComp; }

// Method: getComp
// Usage: p = getComp();
// ---------------------
// Gets the number of each individual component.

int* Ensemble::getComp() { return comp; }

// Method: getVolume
// Usage: getVolume();
// -------------------
// Gets the volume of the simulation box.

double Ensemble::getVolume() { return boxVol; }

// Method: getLength
// Usage: getLength();
// -------------------
// Gets the length of the simulation box.

double Ensemble::getLength() { return boxLen; }

// Method: getPotEnergy
// Usage: getPotEnergy();
// ----------------------
// Gets potential energy following force calculation.

double Ensemble::getPotEnergy() { return potEnergy; }

// Method: getVirial
// Usage: getVirial();
// -------------------
// Gets virial following force calculation.

double Ensemble::getVirial() { return virial; }

// Method: getKineticE
// Usage n = getKineticE();
// ------------------------
// Gets the kinetic energy.

double Ensemble::getKineticE() { return kineticE; }

// Method: getAtoms
// Usage: n = getAtoms();
// ----------------------
// Return array of atoms

Atom** Ensemble::getAtoms() { return atoms; }

// Method:  initAcceleration
// Usage:   initAcceleration();
// ---------------------------
// Performs initialisation on the accelerations of all the atoms
// in the ensemble prior to starting the simulation.
// Requires the force values to have been initialised.

void Ensemble::initAcceleration() {
    int k;
    double *accel, mass;

    accel = new double[3];

    for (int i = 0; i < numAtom; i++) {
        k = 3 * i;
        mass = atoms[i]->getMass();

        for (int j = 0; j < 3; j++)
            accel[j] = atoms[0]->f[k + j] / mass;

        atoms[i]->setAcceleration(&accel[0]);
    }
}

// Method: Ensemble
// Usage: Ensemble;
// ----------------
// Constructor for the Ensemble class

Ensemble::Ensemble(Atom** theAtoms, int nComp, int nAtoms, double temp, double dens, double* mol,
                   int* cmp) {
    atoms = theAtoms;
    numComp = nComp;
    numAtom = nAtoms;
    temperature = temp;
    density = dens;
    molFract = mol;
    comp = cmp;
    potEnergy = 0.0;
    virial = 0.0;
}

// Destructor for Ensemble class
Ensemble::~Ensemble() {}
