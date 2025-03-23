// File:  verletMD.h
// -----------------
// This file contains the definition of the Verlet class,
// which defines all methods and data for the
// Verlet implementation of the integrator method.

#ifndef _verletMD_h_
#define _verletMD_h_

#include "integrtMD.h"

class Verlet : public Integrator {
  public:
    Verlet(Atom**);                                // constructor
    virtual void velVerletP1(int, double, double); // velocity Verlet integrator
    virtual double velVerletP2(int, double, double);
};

#endif
