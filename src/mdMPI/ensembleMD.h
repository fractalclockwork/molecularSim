// File:  ensembleMD.h
// -----------------
// This file contains the definition of abstract class
// Ensemble, which defines all methods and data common to
// different ensemble types (eg derived classes), as well as
// foreshadowing methods of those derived classes through
// the use of abstract methods

#ifndef _ensembleMD_h
#define _ensembleMD_h
#include "atomMD.h"
#include "forceMD.h"

// Class: Ensemble
// ---------------
// The Ensemble class contains the physical attributes
// of the abstract ensemble, from which instantiable ensemble
// classes are derived, such as the NVE ensemble

class Ensemble {
  protected:
    Atom** atoms;       // array of atoms;
    Force* theForce;    // the force
    int numComp;        // number of components
    int numAtom;        // number of atom
    int* comp;          // number of components of each type
    double temperature; // ensemble temperature
    double density;     // ensemble density
    double boxVol;      // ensemble box volume
    double boxLen;      // ensemble box length
    double* molFract;   // mole fraction for each molecule type
    double kineticE;    // Kinetic energy for ensemble
    double potEnergy;   // potential energy of the ensemble
    double virial;      // the virial of the ensemble
  public:
    Ensemble(Atom**, int, int, double, double, double*, int*); // ensemble constructor
    virtual ~Ensemble();                                       // destructor
    Atom** getAtoms();                                         // return the array of atoms
    void setVolume();                                          // determine box volume
    void setLength();                                          // determine box length
    void setComp();          // determine the number of atoms of each type
    int getNumAtom();        // get the total number of atoms
    int* getComp();          // get the number of individual components
    int getNumComp();        // get the total number of components
    double getTemp();        // get the temperature
    double getVolume();      // get the volume of the box
    double getLength();      // get the length of the box
    double getPotEnergy();   // get the potential energy of ensemble
    double getVirial();      // get the virial of the ensemble
    double getKineticE();    // get the kinetic energy of the ensemble
    void initAcceleration(); // initialise the acceleration before loop

    // Abstract Methods
    // ----------------
    // The following methods are defined in the
    // NVEensemble class (nve.h and nve.cpp)

    virtual void initialVelocity() = 0;   // assign the initial velocities
    virtual void initialCoord() = 0;      // place atoms on lattice
    virtual void scaleVelocity() = 0;     // scale velocities to temperature
    virtual void setForce() = 0;          // instigate force calculations
    virtual void setKineticE() = 0;       // set the Kinetic Energy
    virtual void setKineticE(double) = 0; // set the Kinetic Energy;
    virtual double getEnergyLRC() = 0;    // get energy long range correction
    virtual double getVirialLRC() = 0;    // get virial long range correction
    virtual void lrc() = 0;               // long range corrections
};

#endif
