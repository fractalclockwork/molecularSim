// File:  simMC.cpp
// -----------------
// This file contains the implementation of the Simulation class,
// which defines all methods and data for a Molecular Simulation,
// using approaches such as Molecular Dynamics, Monte Carlo etc.

// This code was specifically developed to illustrate concepts in the accompanying book:
// R. J. Sadus, "Molecular Simulation of Fluids: Theory, Algorithms, Object-Orientation,
// and Parallel Computing," 2nd Ed. (Elsevier, Amsterdam, 2023). It can be used freely for
// any not for profit purpose or academic research application. The code has been validated,
// but it would be nonetheless prudent to test it further before publishing any results.
// Check the book's website for any subsequent updates.

#include "simMC.h"
#include <fstream>
#include <iostream>

using namespace std;

// Method:  readSimParameters
// Usage:   readSimParameters();
// -----------------------------
// Reads in parameters for NVE ensemble

void Simulation::readSimParameters() {
    int simulation;
    // open the simMC.dat file
    /*
      ifstream in;
      in.open("simMC.dat");

      if(in.fail()){
             cout << "Cannot open simMC.dat!\n";
             return;
      }
    */

    // Use the DATA_PATH macro defined in CMake
    string dataPath = string(DATA_PATH) + "/simMC.dat";

    cout << "Attempting to open: " << dataPath << endl;

    ifstream in(dataPath);
    if (in.fail()) {
        cout << "Cannot open " << dataPath << "!\n";
        return;
    }

    in >> simulation;
    in.close();

    switch (simulation) {
    case 2: // Monte Carlo
        mc = new MonteCarlo();
        mc->run();
        break;
    default:
        cout << "Invalid simulation selected!  Aborting" << endl;
        break;
    }
    return;
}
