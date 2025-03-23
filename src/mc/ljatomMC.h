// File:  ljatomMC.h
// ----------------
// This file contains the definition of the LJatom class,
// which defines all methods and data for an LJatom
// interacting via the Lennard-Jones 12-6 intermolecular potential.

#ifndef _ljatomMC_h
#define _ljatomMC_h

#include "atomMC.h"

class LJatom : public Atom
{
  private:
    double **epsilon;                  // LJ epsilon for atom pairs
    double **sigma;                    // LJ sigma for atom pairs
  public:
    LJatom(int, double, double **, 
	  double **, double **, int);  // LJ atom constructor
    LJatom();                          // default constructor for array
    virtual void setSigma(double **);  // set sigma for atom pairs
    virtual void setEpsilon(double **);// set epsilon for atom pairs
    virtual double **getEpsilon();     // get epsilon for atom pairs
    virtual double **getSigma();       // get sigma for atom pairs
};

#endif
