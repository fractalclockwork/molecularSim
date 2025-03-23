// File:  ensembleMC.h
//  -----------------
// This file contains the definition of abstract class
// Ensemble, which defines all methods and data common to
// different ensemble types (eg derived classes), as well as
// foreshadowing methods of those derived classes through
// the use of abstract methods

#ifndef _ensembleMC_h
#define _ensembleMC_h

#include "atomMC.h"
#include "energyMC.h"

// Class: Ensemble
// ---------------
// The Ensemble class contains the physical attributes
// of the abstract ensemble, from which instantiable ensemble
// classes are derived, such as the NVE ensemble

class EnsembleMC {
  protected:
    AtomMC** atoms;    // array of atoms
    Energy* theEnergy; // reference to energy
    int numComp;       // number of components
    int numAtom;       // number of atom
    int* comp;         // number of components of each type
    double totEfixed;  // total fixed energy
    double density;    // ensemble density
    double boxVol;     // ensemble box volume
    double boxLen;     // ensemble box length
    double* molFract;  // mole fraction for each molecule type
    double kineticE;   // Kinetic energy for ensemble
    double potEnergy;  // potential energy of the ensemble
    double trialPotE;

  public:
    EnsembleMC(AtomMC**, int, int, double, double, double*, int*); // ensemble constructor
    virtual ~EnsembleMC();                                         // destructor
    AtomMC** getAtoms();                                           // return the array of atoms
    void setVolume();                                              // determine box volume
    void setLength();                                              // determine box length
    void setComp(); // determine the number of atoms of each type
    void setKineticE(double);
    int getNumAtom();      // get the total number of atoms
    int* getComp();        // get the number of individual components
    int getNumComp();      // get the total number of components
    double getVolume();    // get the volume of the box
    double getLength();    // get the length of the box
    double getPotEnergy(); // get the potential energy of ensemble
    double getKineticE();  // get the kinetic energy of the ensemble
    double gettotEfixed(); // get the total fixed energy of the ensemble

    virtual double getTrialPotE(int, AtomMC**) = 0;
    virtual void reSetEnergy(double) = 0;
    virtual void initialCoord() = 0;   // place atoms on lattice
    virtual void setEnergy() = 0;      // Order N*N energy calculation
    virtual double getEnergyLRC() = 0; // get energy long range correction
    virtual void lrc() = 0;            // trigger the long range corrections
};

#endif
