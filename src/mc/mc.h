// File:  mc.h
// -----------------
// This file contains the definition of the MonteCarlo class,
// which defines all methods and data for a Monte Carlo approach
// to molecular simulation.

#ifndef _mc_h_
#define _mc_h_

#include "ensembleMC.h"
#include "ljatomMC.h"

class MonteCarlo {
  private:
    Atom** atom;
    Ensemble* ensemble;
    int theEnsemble;

  protected:
    int nSize, nEquil, nStep;

  public:
    MonteCarlo();
    int getnSize();
    int getnEquil();
    int getnStep();
    void runNVE();
    void readInNVE(int);
    void run();
};

#endif
