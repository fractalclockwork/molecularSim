// File:  gpcMD.h
// -----------------
// This file contains the definition of the GearPC class,
// which defines all methods and data for the Gear Predictor
// Corrector implementation of the integrator method.

#ifndef _gpcMD_h_
#define _gpcMD_h_

#include "integrtMD.h"

class GearPC : public Integrator {
  public:
    GearPC(Atom**);                                // constructor
    virtual void gearPredict(int, double, double); // Gear's predictor
    virtual double gearCorrect(int, double, double);
};

#endif
