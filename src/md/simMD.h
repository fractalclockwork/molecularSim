// File:  simMD.h
// -----------------
// This file contains the definition of the Simulation class,
// which defines all methods and data for a Molecular Simulation,
// using approaches such as Molecular Dynamics, Monte Carlo etc.

#ifndef _simMD_h_
#define _simMD_h_

#include "md.h"

class Simulation
{
  private:
    MolecularDynamics *md;
  public:
    void readSimParameters();
    
};

#endif
