// File:  ljenergyMC.h
//  -------------
// This file contains the definition of the LJenergy class,
// which defines all methods and data for a Energy object
// which interacts according to the Lennard-Jones
// intermolecular potential.
 
#ifndef _ljenergyMC_h_
#define _ljenergyMC_h_

#include "energyMC.h"

class LJenergy : public Energy 
{
  private:
   double **epsilon, **sigma, **rCut;  // reference to the epsilon,
				       // sigma, and rCut values
  public:
   LJenergy(AtomMC **);		       // constructor for LJenergy
   virtual void setEnergy(int, double,
            double * );                // N*N order energy calculation
   virtual double getTrialPotE(int, double, int, AtomMC **);	       
   virtual void lrc(int, int *, double,
		 double *);            // long range corrections
};

#endif

