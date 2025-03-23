// File:  ljforceMD.h
//  -------------
// This file contains the definition of the LJforce class,
// which defines all methods and data for a Force object
// which interacts according to the Lennard-Jones
// intermolecular potential.

#ifndef _ljforceMD_h_
#define _ljforceMD_h_

#include "forceMD.h"

class LJforce : public Force {
  private:
    double **epsilon, **sigma, **rCut; // reference to LJ parameters and rCut
  public:
    LJforce(Atom**);                                       // constructor for LJforce
    virtual void setForce(int, double, double*, double*);  // set force calculations
    virtual void lrc(int, int*, double, double*, double*); // long range corrections
};

#endif
