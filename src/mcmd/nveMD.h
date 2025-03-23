// File:  nveMD.h
// -----------------
// This file contains the definition of the NVEensemble class,
// which defines all methods and data for the microcanonical ensemble.

#ifndef _nveMD_h
#define _nveMD_h

#include "ensembleMD.h"

class NVEensemble : public Ensemble {
  private:
    double energyLRC; // long range energy correction
    double virialLRC; // long range virial correction
  public:
    NVEensemble(Atom**, int, int, double, double, double*, int*); // ensemble constructor
    void readInNVE();                                             // read in ensemble data
    virtual void initialVelocity();   // assign initial velocities for the atoms
    virtual void initialCoord();      // place atoms on lattice
    virtual void scaleVelocity();     // scale velocities to temperature
    virtual void setForce();          // instigate force calculations
    virtual void setKineticE();       // set the ensemble kinetice nergy
    virtual void setKineticE(double); // set the ensmble kinetic energy
    virtual double getEnergyLRC();    // get energy long range correction
    virtual double getVirialLRC();    // get virial long range correction
    virtual void lrc();               // trigger the long range corrections
};

#endif
