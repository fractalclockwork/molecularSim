// File:  intergrtMD.h
//  -----------------
// This file contains the definition of the abstract class Integrator,
// which defines all methods and data common to intergrator methods.

#ifndef _integrtMD_h_
#define _integrtMD_h_

#include "atomMD.h"

class Integrator {
  protected:
    Atom** atoms;

  public:
    Integrator(Atom**);                                  // constructor for integrator
    virtual ~Integrator();                               // destructor for integrator
    virtual void gearPredict(int, double, double) = 0;   // Gear's predictor
    virtual double gearCorrect(int, double, double) = 0; // Gear's Corrector
};

#endif
