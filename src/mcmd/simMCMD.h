// File:  simMCMD.h
// -----------------
// This file contains the definition of the Simulation class,
// which defines all methods and data for a Molecular Simulation,
// using approaches such as Molecular Dynamics, Monte Carlo etc.

#ifndef _simMCMD_h_
#define _simMCMD_h_

#include "mc.h"
#include "md.h"

class Simulation {
  private:
    MolecularDynamics* md;
    MonteCarlo* mc;

  public:
    void readSimParameters();
};

#endif
