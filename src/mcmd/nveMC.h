//  File:  nveMC.h
//  -----------------
//  This file contains the definition of the NVEensemble class,
//  which defines all methods and data for the microcanonical ensemble.

#ifndef _nveMC_h
#define _nveMC_h


#include "ensembleMC.h"

class NVEensembleMC : public EnsembleMC
{
  private:
    double energyLRC;    // long range energy correction
  public:
    NVEensembleMC(AtomMC **, int, int, double,
	 double, double *, int *);// ensemble constructor
    void readInNVE();                   // read in ensemble data
    virtual void initialCoord();        //place atoms on lattice

    virtual double getTrialPotE(int, AtomMC **);
    virtual void reSetEnergy(double);
    virtual void setEnergy();           // Order N*N energy calculation
    virtual double getEnergyLRC();      // get energy long range correction
    virtual void lrc();                 // trigger the long range corrections
};

#endif
