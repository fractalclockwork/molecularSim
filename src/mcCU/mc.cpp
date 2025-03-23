// File:  mc.cpp
// Class: MonteCarlo
// -----------------
// This file contains the implementation details for the Monte Carlo
// class to perform an MC simulation

// This code was specifically developed to illustrate concepts in the accompanying book:
// R. J. Sadus, "Molecular Simulation of Fluids: Theory, Algorithms, Object-Orientation,
// and Parallel Computing," 2nd Ed. (Elsevier, Amsterdam, 2023). It can be used freely for
// any not for profit purpose or academic research application. The code has been validated,
// but it would be nonetheless prudent to test it further before publishing any results.
// Check the book's website for any subsequent updates.

#include "mc.h"
#include "nveMC.h"
#include "auxfunc.h"
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
using namespace std;

void MonteCarlo::runNVE()
{
  ofstream out;
  out.open("mc_nve.out", ios::out);

  if(out.fail()) {
     cout << "could not open output file!" << endl;
     return;
  }
  int i, index, j = 0, k, l, n, nTotal, seed = -50, counter = 0, accept = 0;
  int nStep  = getnStep();
  int nEquil = getnEquil();
  int nSize  = getnSize();
  int num    = ensemble->getNumAtom();
  double length = ensemble->getLength();
  double avPotEnergy = 0.0;
  double avKinEnergy = 0.0;
  double avTotEnergy = 0.0;
  double avTemperature = 0.0;
  double erKinEnergy = 0.0;
  double erPotEnergy = 0.0;
  double erTotEnergy = 0.0;
  double erTemperature = 0.0;
  double *accumPotE, *accumKinE;
  double *accumTotE, *accumTemp;
  double *trial, *position, f;
  double **eps, **sig, **cut;
  double eOld, eNew, eTrial, eDelta, trialPotE, potE, kinE, totE, temp;
  double rMax, eFixed, kNew, kOld, wFactor;
  double tally = 0.0, lrc;
  double const MAX = 0.1;
  rMax = MAX*length; 
  nTotal = nStep - nEquil;
  n = nTotal/nSize;
  f = 3.0*num/2.0 - 1.0; 
  accumPotE = new double [n];
  accumKinE = new double [n];
  accumTotE = new double [n];
  accumTemp = new double [n];
  position = new double [3];
  trial = new double [3];

  for (i = 0; i < n; i++)
  {
    accumPotE[i] = 0.0;
    accumKinE[i] = 0.0;
    accumTotE[i] = 0.0;
    accumTemp[i] = 0.0;
  }
  
  atom = ensemble->getAtoms();
  ensemble->lrc(); //perform long range corrections
  lrc = ensemble->getEnergyLRC();

  // Initialization of 'master'arrays for GPUs, which
  // requires defining public arrays in the atom class. The alterative
  // of protected variables accessed via methods is computationally
  // costly when using GPUs with the current code design..
  
  for(i = 0; i < num; i++)
  {
   position = atom[i]->getPosition();
   // A single vector array is only used to avoid separate arrays
   // for the x, y and z components 
   atom[0]->r[3*i]     = position[0];
   atom[0]->r[3*i + 1] = position[1];
   atom[0]->r[3*i + 2] = position[2];

   l = atom[0]->kind[i] = atom[i]->getType();

   eps = atom[0]->getEpsilon();
   sig = atom[0]->getSigma();
   cut = atom[0]->getrCutOff();
   
   // The transformation  below is used to avoid passing 2-D 
   // arrays to the kernels. The unlike contributions for mixtures are
   // calculated in the kernels using combining rules Aadditional
   // i-j arrays could be used as an alternative to combining rules.
   atom[0]->epsii[i]  = eps[l][l];
   atom[0]->sigii[i]  = sig[l][l];
   atom[0]->rCutii[i] = cut[l][l];

  }

  ensemble->setEnergy();  //set initial energy
  potE = ensemble->getPotEnergy();
  eFixed = ensemble->gettotEfixed();
  kNew   = eFixed - (potE + lrc); //determine inital kinetic energy

  if(kNew > 0)
      ensemble->setKineticE(kNew);
  else
  {
      cout << "Initial configutaion has negative kinetic energy. Aborting." << endl;
  }
  cout << "Results will be directed to file \"mc_nve.out\"" << endl;

  for(i = 0; i < nStep; i++)
  {
    
     for (k = 0; k < num; ++k)
     {
       ++counter;
       potE = ensemble->getPotEnergy();
       kOld = ensemble->getKineticE();
         
       //Select atom randomly
       index = (int) num*random(&seed);
       position  = atom[index]->getPosition();
       atom[index]->setTrialPosition(position);
       
       // determine energy at current position
       eOld = ensemble->getTrialPotE(index,atom);

       // find trial position
       trial[0] =  position[0] + rMax*(2.0*random(&seed)-1.0);
       trial[1] =  position[1] + rMax*(2.0*random(&seed)-1.0);
       trial[2] =  position[2] + rMax*(2.0*random(&seed)-1.0);

       // apply periodic boundary conditions
       trial[0] -= length * nearestInt(trial[0],length);
       trial[1] -= length * nearestInt(trial[1],length);
       trial[2] -= length * nearestInt(trial[2],length);
       atom[index]->setTrialPosition(trial);
       
       // calculate energy at trial position
       eNew = ensemble->getTrialPotE(index,atom);
       eDelta = eNew - eOld;

       trialPotE = potE + eDelta;
       kNew = eFixed - (trialPotE +lrc);
  
       if (kNew > 0)
       {
         wFactor = pow(kNew/kOld, f);
         if (wFactor > 1 || wFactor > random(&seed))  //accept move
         {
           atom[index]->reSetPosition();

           // Update the 'master' position vector array used for GPUs
           atom[0]->r[3*index]      = trial[0];
           atom[0]->r[3*index + 1]  = trial[1];
           atom[0]->r[3*index + 2]  = trial[2];

           ensemble->setKineticE(kNew);
           ensemble->reSetEnergy(trialPotE);
           accept++;
 
         }
        }
      }
      
      tally = accept/(double) counter;
       if (i < nEquil)
       {
          if (tally < 0.5)
            rMax *=0.95;
          else
            rMax *= 1.05;
        }
      
      if(i >= nEquil)
      {
         potE = ensemble->getPotEnergy() + lrc;
         kinE = ensemble->getKineticE();
         totE = potE + kinE;
         temp = 2.0*kinE/(3.0*num);
         avPotEnergy += potE;
         avKinEnergy += kinE;
         avTotEnergy += totE;
         avTemperature += temp;

        if((i != nEquil) && ((i % nSize) == 0))
        {
          accumPotE[j] /= nSize;
	      accumKinE[j] /= nSize;
          accumTotE[j] /= nSize;
          accumTemp[j] /= nSize;
	      j++;
	    }

         accumPotE[j] += potE;
         accumKinE[j] += kinE;
         accumTotE[j] += totE;
         accumTemp[j] += temp;
      }
  }
 
  avPotEnergy /= nTotal;
  avKinEnergy /= nTotal;
  avTotEnergy /= nTotal;
  avTemperature /= nTotal;

  for(i = 0; i < n; i++)
  {
    erPotEnergy += pow((accumPotE[i] - avPotEnergy),2);
    erKinEnergy += pow((accumKinE[i] - avKinEnergy),2);
    erTotEnergy += pow((accumTotE[i] - avTotEnergy),2);
    erTemperature += pow((accumTemp[i] - avTemperature),2); 
  }

  out <<"Average Potential Energy:/N\t" << avPotEnergy/num
       <<"  +/-  "<<sqrt(erPotEnergy)/(nTotal*num)<<endl;
  out <<"Average Kinetic Energy/N:\t" << avKinEnergy/num
       <<"  +/-  "<<sqrt(erKinEnergy)/(nTotal*num)<<endl;
  out <<"Average Total Energy/N:\t\t" << avTotEnergy/num
      <<"  +/-  "<<sqrt(erTotEnergy)/(nTotal*num)<<endl;
  out <<"Average Temperature:\t\t" << avTemperature
      <<"  +/-  "<<sqrt(erTemperature)/nTotal<<endl; 
  out <<"Acceptance Rate:\t\t"<<100*tally <<" %"<<endl;
  out.close();
}

// Method: readInNVE
// Usage: readInNVE();
// ------------------- 
// ReadIn reads-in the NVE ensemble settings from the
// NVEfileMC.dat data file.

void MonteCarlo::readInNVE(int interPotential)
{
  int i, j;
  double numType, *molFract, **epsilon, **sigma, 
         **rCut,  potE, density;
  int dimensions = 3, numComp, numAtom, *comp;
  Atom *newAtom, **atoms;

  // open the NVEfileMC.dat file
/*
  ifstream in;
  in.open("NVEfileMC.dat");

  if(in.fail()){
	 cout << "Cannot open NVEFileMC.dat!\n";
	 return;
  }
*/
    // Use the DATA_PATH macro defined in CMake
    string dataPath = string(DATA_PATH) + "/NVEfileMC.dat";

    cout << "Attempting to open: " << dataPath << endl;

    ifstream in(dataPath);
    if (in.fail()) {
        cout << "Cannot open " << dataPath << "!\n";
        return;
    } 

  in >> numComp;

  if(numComp <= 0){
     cout << "Number of components must be > 0!\n";
	 return;
  }

  comp = new int[numComp];

  in >> numAtom;

  if(numAtom < 4){
    cout << "At least 4 atoms are required!\n";
    return;
  }

  if(!(molFract = new double [numComp]))
  {
    cout << "Cannot allocate memory to molFract" << endl;
	 return;
  }
  for(i = 0; i < numComp; i++)
      in >> molFract[i];

  double *mass = new double[numComp];

  for(i = 0; i < numComp; i++)
	 in >> mass[i];

  in >> potE;
  in >> density;
  potE *= numAtom;
  in.close();

  // construct array of atoms
  if(interPotential == 1) // Lennard-Jones
  {
/*
	in.open("paramLJ.dat");
	if(in.fail()){
	  cout << "Cannot open paramLJ.dat!\n";
	  return;
	}
*/
    // Use the DATA_PATH macro defined in CMake
    string dataPath = string(DATA_PATH) + "/paramLJ.dat";

    cout << "Attempting to open: " << dataPath << endl;

    ifstream in(dataPath);
    if (in.fail()) {
        cout << "Cannot open " << dataPath << "!\n";
        return;
    } 

	// assign memory to arrays
	if(!(atoms = new Atom * [numAtom]))
	{
	 cout << "Cannot allocate memory to atoms!\n";
	 return;
	}

	epsilon = getMemory(numComp);
	sigma = getMemory(numComp);
	rCut = getMemory(numComp);

	for (i = 0; i < numComp; i++)
	{
	  for (j = 0; j < numComp; j++)
	  {
	    in >> epsilon[i][j] >> sigma[i][j] >> rCut[i][j];

	    //adjust for indistinguishable pairs
	    if (j!=i)
	    {
	     epsilon[j][i] = epsilon[i][j];
	     sigma[j][i] = sigma[i][j];
	     rCut[j][i] = rCut[i][j];
	    }
	  }
	}

	// close input file
	in.close();

	int *atnum = new int[numComp + 1];
	atnum[0] = 0;

	// calculate the numebr of atoms of each type from the molFract
	for(i = 0; i < numComp; i++)
	{
	  numType = numAtom * molFract[i];
	  int rem = (int) numType;

	  if ((numType - rem) >= .5)
	     atnum[i+1] = (int) numType + 1;
	  else
	     atnum[i+1] = (int) numType;
	}

	int sum = 0;
	for (i = 0; i <= numComp; i++)
	sum += atnum[i];

	while (atnum[numComp] < numAtom)
	       atnum[numComp]++;

	for(i = 0; i < numComp; i++)
		atnum[i+1] += atnum[i];

	// loop through the number of components
	for (i = 0; i < numComp; i++)
	{
	  //loop through the number of atoms of that type
	  for (int j = atnum[i]; j < atnum[i+1]; j++)
	  {
	    //create atom(s) and assign mass, type, sigma, epsilon, and rCut
	    //(derivatives - 2) because acceleration and velocity are stored
	    //separately from the higher derivatives of time
	    newAtom = new LJatom(i, mass[i], epsilon, 
				 sigma, rCut, dimensions, numAtom);

	    //store reference to atom in array
	    atoms[j] = newAtom;
  	  }
	}
  } // end of Lennard-Jones atom array creation

  ensemble = new NVEensemble(atoms, numComp, numAtom, potE,
                             density, molFract, comp);

  return;
}

// Method: getnSize()
// Usage: n = getSize();
// --------------------
// Return the size of the blocks for averaging
// ensemble properties.
 	
int MonteCarlo::getnSize()
{
   return nSize;
}

// Method: getnEquil()
// Usage: n = getnEquil();
// -----------------------
// Return the number of steps prior to
// equilibrium.

int MonteCarlo::getnEquil()
{
  return nEquil;
}

// Method: getnStep()
// Usage: n = getnStep();
// ----------------------
// Return the total number of simulation
// steps.

int MonteCarlo::getnStep()
{
   return nStep;
}
 
// Method:  run
// Usage:  run();
// --------------
// Runs the appropriate simularion run method 
// determined by the choice of ensemble.

void MonteCarlo::run()
{
  switch(theEnsemble)
  {
    case 1: // NVE ensemble
      runNVE();
      break;
      // other ensembles can be inserted here
  }
}

// Constructor
// -----------
// Accesses parameter file mc.dat, which identifies the choice of
// the choice of intermolecular potential, and ensemble, and
// constructs the MonteCarlo object.
 
MonteCarlo::MonteCarlo()
{
  int potential;
/*
  ifstream in;

  in.open("mc.dat");

  if(in.fail())
  {
    cout << "Unable to open mc.dat for MC parameters" << endl;
    return;
  }
*/
    // Use the DATA_PATH macro defined in CMake
    string dataPath = string(DATA_PATH) + "/mc.dat";

    cout << "Attempting to open: " << dataPath << endl;

    ifstream in(dataPath);
    if (in.fail()) {
        cout << "Cannot open " << dataPath << "!\n";
        return;
    } 

  in >> nStep;
  in >> nEquil;
  in >> nSize;
  in >> theEnsemble >> potential;
  in.close();

  switch(theEnsemble)
  {
    case 1: //NVE ensemble
      readInNVE(potential);
      break;
     //other ensembles here
  }

}

