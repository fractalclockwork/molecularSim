// File:  forceMD.h
// -----------------
// This file contains the definition of abstract class
// Force, which defines all methods and data common to
// different force types (eg derived classes), as well as
// foreshadowing methods of those derived classes through
// the use of abstract methods

#ifndef _forceMD_h
#define _forceMD_h

#include "atomMD.h"

class Force {
  protected:
    Atom** atoms; // reference to the array of atoms
  public:
    Force(Atom**);    // constructor for Force
    virtual ~Force(); // destructor for Force

    // Abstract Methods
    // ----------------
    // abstract methods for subclasses
    // for LJforce
    virtual void setForce(int, double, double*, double*) = 0;  // set force calculations
    virtual void lrc(int, int*, double, double*, double*) = 0; // long range corrections
};

#endif
