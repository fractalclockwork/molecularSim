//  File:  simMC.h
//  -----------------
//  This file contains the definition of the Simulation class,
//  which defines all methods and data for a Molecular Simulation,
//  using approaches such as Molecular Dynamics, Monte Carlo etc.

#ifndef _simMC_h_
#define _simMC_h_

#include "mc.h"

class Simulation{
  private:
     MonteCarlo *mc;
  public:
    void readSimParameters();
};

#endif
