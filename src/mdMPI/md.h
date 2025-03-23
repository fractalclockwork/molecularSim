// File:  md.h
// -----------------
// This file contains the definition of the MolecularDynamics class,
// which defines all methods and data for a Molecular Dynamics approach
// to molecular simulation

#ifndef _md_h_
#define _md_h_

#include "ensembleMD.h"
#include "integrtMD.h"

class MolecularDynamics {
  private:
    Atom** atom;
    Ensemble* ensemble;
    Integrator* integrator;
    int theEnsemble, derivatives;

  protected:
    int nSize, nEquil, nStep;
    double tStep;

  public:
    MolecularDynamics();
    int getnSize();
    int getnEquil();
    int getnStep();
    double gettStep();
    void runNVE();
    void readInNVE(int);
    void run();
};

#endif
