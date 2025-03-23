// File:  md.cpp
// Class: MolecularDynamics
// ------------------------
// This file contains the implementation details for the Molecular dynamics
// class to perform an MD simulation.

// This code was specifically developed to illustrate concepts in the accompanying book:
// R. J. Sadus, "Molecular Simulation of Fluids: Theory, Algorithms, Object-Orientation,
// and Parallel Computing," 2nd Ed. (Elsevier, Amsterdam, 2023). It can be used freely for
// any not for profit purpose or academic research application. The code has been validated,
// but it would be nonetheless prudent to test it further before publishing any results.
// Check the book's website for any subsequent updates.

#include "md.h"
#include "auxfunc.h"
#include "gpcMD.h"
#include "ljatomMD.h"
#include "nveMD.h"
#include "verletMD.h"
#include <fstream>
#include <iostream>
#include <math.h>
#include <mpi.h>

using namespace std;

void MolecularDynamics::runNVE() {
    ofstream out;
    out.open("md_nve.out", ios::out);

    if (out.fail()) {
        cout << "could not open output file!" << endl;
        return;
    }

    int i, j = 0, k, n, nTotal;
    int nStep = getnStep();
    int nEquil = getnEquil();
    int nSize = getnSize();
    int num = ensemble->getNumAtom();
    double length = ensemble->getLength();
    double volume = ensemble->getVolume();
    double tStep = gettStep(), KinE;
    double potE, virial, kinE, totE, temp, lrcE, lrcW;
    double avPotEnergy = 0.0;
    double avKinEnergy = 0.0;
    double avTotEnergy = 0.0;
    double avVirial = 0.0;
    double avTemperature = 0.0;
    double erKinEnergy = 0.0;
    double erPotEnergy = 0.0;
    double erTotEnergy = 0.0;
    double erVirial = 0.0;
    double erTemperature = 0.0;
    double pressure, erPress;
    double* position;
    double *accumPotE, *accumKinE;
    double *accumTotE, *accumVir, *accumTemp;

    nTotal = nStep - nEquil;
    n = nTotal / nSize;
    accumPotE = new double[n];
    accumKinE = new double[n];
    accumTotE = new double[n];
    accumVir = new double[n];
    accumTemp = new double[n];

    for (i = 0; i < n; i++) {
        accumPotE[i] = 0.0;
        accumKinE[i] = 0.0;
        accumTotE[i] = 0.0;
        accumVir[i] = 0.0;
        ;
        accumTemp[i] = 0.0;
    }

    ensemble->lrc();                 // perform long range corrections
    lrcE = ensemble->getEnergyLRC(); // long range correction for energy
    lrcW = ensemble->getVirialLRC(); // long range correction for virial

    // Iniialise 'master' root node position array to facilitate MPI implementaion

    atom = ensemble->getAtoms();

    for (i = 0; i < num; i++) {
        position = atom[i]->getPosition();
        k = 3 * i;
        atom[0]->r[k] = position[0];
        atom[0]->r[k + 1] = position[1];
        atom[0]->r[k + 2] = position[2];
    }

    ensemble->setForce();         // set initial force prior to setting acceleration
    ensemble->initAcceleration(); // initialise acceleration
    ensemble->initialVelocity();  // initialise velocity
    ensemble->setKineticE();      // determine kinetic energy
    ensemble->scaleVelocity();    // scale initial velocities
    ensemble->setKineticE();      // re-determine kinetic energy with scaled velcities

    cout << "Results will be directed to the file \"md_nve.out\"" << endl;

    for (i = 0; i < nStep; i++) {
        // Predict new position etc, then determine the new force and KE to update
        // everthing in the corrector step

        integrator->gearPredict(num, length, tStep);
        ensemble->setForce();
        kinE = integrator->gearCorrect(num, length, tStep);
        ensemble->setKineticE(kinE);

        // Periodically scale velocities to obtain the desired temperature

        if (((i % 10) == 0) && i < nEquil)
            ensemble->scaleVelocity();

        if (i >= nEquil) {
            potE = ensemble->getPotEnergy() + lrcE;
            kinE = ensemble->getKineticE();
            virial = ensemble->getVirial() + lrcW;
            totE = potE + kinE;
            temp = 2 * kinE / (3 * num - 3);
            avPotEnergy += potE;
            avKinEnergy += kinE;
            avTotEnergy += totE;
            avVirial += virial;
            avTemperature += temp;

            if ((i != nEquil) && ((i % nSize) == 0)) {
                accumPotE[j] /= nSize;
                accumKinE[j] /= nSize;
                accumTotE[j] /= nSize;
                accumVir[j] /= nSize;
                accumTemp[j] /= nSize;
                j++;
            }
            accumPotE[j] += potE;
            accumKinE[j] += kinE;
            accumTotE[j] += totE;
            accumVir[j] += virial;
            accumTemp[j] += temp;
        }
    }

    avPotEnergy /= nTotal;
    avKinEnergy /= nTotal;
    avTotEnergy /= nTotal;
    avVirial /= nTotal;
    ;
    avTemperature /= nTotal;

    for (i = 0; i < n; i++) {
        erPotEnergy += pow((accumPotE[i] - avPotEnergy), 2);
        erKinEnergy += pow((accumKinE[i] - avKinEnergy), 2);
        erTotEnergy += pow((accumTotE[i] - avTotEnergy), 2);
        erVirial += pow((accumVir[i] - avVirial), 2);
        ;
        erTemperature += pow((accumTemp[i] - avTemperature), 2);
    }

    pressure = avVirial / volume + num * ensemble->getTemp() / volume;
    erPress = erVirial / volume;

    out << "Potential Energy/N:\t" << avPotEnergy / num << "  +/-  "
        << sqrt(erPotEnergy) / (nTotal * num) << endl;
    out << "Kinetic Energy/N:\t " << avKinEnergy / num << "  +/-  "
        << sqrt(erKinEnergy) / (nTotal * num) << endl;
    out << "Total Energy/N:\t\t" << avTotEnergy / num << "  +/-  "
        << sqrt(erTotEnergy) / (nTotal * num) << endl;
    out << "Pressure:\t\t " << pressure << "  +/-  " << sqrt(erPress) / (nTotal * num) << endl;
    out << "Temperature:\t\t " << avTemperature << "  +/-  " << sqrt(erTemperature) / nTotal
        << endl;
    out.close();
}

// Method: readInNVE
// Usage: readInNVE();
// ---------------------
// ReadIn reads-in the NVE ensemble settings from the
// NVEfileMD.dat data file.

void MolecularDynamics::readInNVE(int interPotential) {
    int i, j;
    double numType, *molFract, **epsilon, **sigma, **rCut;
    double temperature, density;
    int dimensions = 3, numComp, numAtom, *comp;
    Atom *newAtom, **atoms;

    // open the NVEfileMD.dat file
    /*
      ifstream in;
      in.open("NVEfileMD.dat");

      if(in.fail())
      {
       cout << "Cannot open NVEFileMD.dat!\n";
       return;
      }
    */
    // Use the DATA_PATH macro defined in CMake
    string dataPath = string(DATA_PATH) + "/NVEfileMD.dat";

    cout << "Attempting to open: " << dataPath << endl;

    ifstream in(dataPath);
    if (in.fail()) {
        cout << "Cannot open " << dataPath << "!\n";
        return;
    }

    in >> numComp;

    if (numComp <= 0) {
        cout << "Number of components must be > 0!\n";
        return;
    }

    comp = new int[numComp];

    in >> numAtom;

    if (numAtom < 4) {
        cout << "At least 4 atoms are required!\n";
        return;
    }

    if (!(molFract = new double[numComp])) {
        cout << "Cannot allocate memory to molFract" << endl;
        return;
    }

    for (i = 0; i < numComp; i++)
        in >> molFract[i];

    double* mass = new double[numComp];

    for (i = 0; i < numComp; i++)
        in >> mass[i];

    in >> temperature;
    in >> density;
    in.close();

    // construct array of atoms
    if (interPotential == 1) // Lennard-Jones
    {
        /*
            in.open("paramLJ.dat");
            if(in.fail())
            {
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
        if (!(atoms = new Atom*[numAtom])) {
            cout << "Cannot allocate memory to atoms!\n";
            return;
        }

        epsilon = getMemory(numComp);
        sigma = getMemory(numComp);
        rCut = getMemory(numComp);

        for (i = 0; i < numComp; i++) {
            for (j = 0; j < numComp; j++) {
                in >> epsilon[i][j] >> sigma[i][j] >> rCut[i][j];

                // adjust for indistinguishable pairs
                if (j != i) {
                    epsilon[j][i] = epsilon[i][j];
                    sigma[j][i] = sigma[i][j];
                    rCut[j][i] = rCut[i][j];
                }
            }
        }

        // close input file
        in.close();

        int* atnum = new int[numComp + 1];
        atnum[0] = 0;

        // calculate the numebr of atoms of each type from the molFract
        for (i = 0; i < numComp; i++) {
            numType = numAtom * molFract[i];
            int rem = (int)numType;
            if ((numType - rem) >= .5)
                atnum[i + 1] = (int)numType + 1;
            else
                atnum[i + 1] = (int)numType;
        }

        int sum = 0;
        for (i = 0; i <= numComp; i++)
            sum += atnum[i];

        while (atnum[numComp] < numAtom)
            atnum[numComp]++;

        for (i = 0; i < numComp; i++)
            atnum[i + 1] += atnum[i];

        // loop through the number of components
        for (i = 0; i < numComp; i++) {
            // loop through the number of atoms of that type
            for (int j = atnum[i]; j < atnum[i + 1]; j++) {
                // create atom(s) and assign mass, type, sigma, epsilon, and rCut
                //(derivatives - 2) because acceleration and velocity are stored
                newAtom = new LJatom(i, mass[i], epsilon, sigma, rCut, dimensions,
                                     (derivatives - 2), numAtom);

                // store reference to atom in array
                atoms[j] = newAtom;
            }
        }
    } // end of Lennard-Jones atom array creation

    ensemble = new NVEensemble(atoms, numComp, numAtom, temperature, density, molFract, comp);

    return;
}

// Method: getnSize;
// Usage:  n = getSize();
// ----------------------
// Return the size of the blocks for gathering statistics

int MolecularDynamics::getnSize() { return nSize; }

// Method: getnEquil
// Usage:  n = getnEquil();
// ------------------------
// Return the equilibration period

int MolecularDynamics::getnEquil() { return nEquil; }

// Method: getnStep
// Usage: n = getnStep();
// -----------------------
// Return the number of time steps or cycles

int MolecularDynamics::getnStep() { return nStep; }

// Method: gettStep
// Usage: n = gettStep();
// ----------------------
// Return the time step

double MolecularDynamics::gettStep() { return tStep; }

// Method:  run
// Usage:  run();
// --------------
// Runs the appropriate simularion run method
// determined by the choice of ensemble.

void MolecularDynamics::run() {
    switch (theEnsemble) {
    case 1: // NVE ensemble
        runNVE();
        break;
        // other ensembles can be inserted here
    }
}

// Constructor
// -----------
// Accesses parameter file md.dat, which identifies the choice of
// intermolecular potential, ensemble, and integrator method, and
// constructs the MolecularDynamics object.

MolecularDynamics::MolecularDynamics() {
    int integratorM, interPotential;
    /*
      ifstream in;

      in.open("md.dat");

      if(in.fail())
      {
        cout << "Unable to open md.dat for MD parameters" << endl;
        return;
      }
    */
    // Use the DATA_PATH macro defined in CMake
    string dataPath = string(DATA_PATH) + "/md.dat";

    cout << "Attempting to open: " << dataPath << endl;

    ifstream in(dataPath);
    if (in.fail()) {
        cout << "Cannot open " << dataPath << "!\n";
        return;
    }

    in >> nStep;
    in >> nEquil;
    in >> nSize;
    in >> tStep;

    in >> integratorM;
    if (integratorM == 1) // gearPredictorCorrector
        in >> derivatives;

    in >> theEnsemble >> interPotential;
    in.close();

    switch (theEnsemble) {
    case 1: // NVE ensemble
        readInNVE(interPotential);
        break;
        // other ensembles here
    }

    switch (integratorM) {
    case 1:
        integrator = new GearPC(ensemble->getAtoms());
        //      integrator = new Verlet(ensemble->getAtoms());
        break;
        // other intergrator methods can be inserted here
    }
}
