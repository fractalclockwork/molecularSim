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
#include "auxfunc.h"
#include "nveMC.h"
#include <fstream>
#include <iostream>
#include <math.h>

#include <toml++/toml.h>

using namespace std;

void MonteCarlo::runNVE() {
    ofstream out;
    out.open("mc_nve.out", ios::out);

    if (out.fail()) {
        cout << "could not open output file!" << endl;
        return;
    }
    int i, index, j = 0, k, n, nTotal, seed = -50, counter = 0, accept = 0;
    int nStep = getnStep();
    int nEquil = getnEquil();
    int nSize = getnSize();
    int num = ensemble->getNumAtom();
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
    double eOld, eNew, eTrial, eDelta, trialPotE, potE, kinE, totE, temp;
    double rMax, eFixed, kNew, kOld, wFactor;
    double tally = 0.0, lrc;
    double const MAX = 0.1;
    rMax = MAX * length;
    nTotal = nStep - nEquil;
    n = nTotal / nSize;
    f = 3.0 * num / 2.0 - 1.0;
    accumPotE = new double[n];
    accumKinE = new double[n];
    accumTotE = new double[n];
    accumTemp = new double[n];
    position = new double[3];
    trial = new double[3];

    for (i = 0; i < n; i++) {
        accumPotE[i] = 0.0;
        accumKinE[i] = 0.0;
        accumTotE[i] = 0.0;
        accumTemp[i] = 0.0;
    }

    atom = ensemble->getAtoms();
    ensemble->lrc(); // perform long range corrections
    lrc = ensemble->getEnergyLRC();
    ensemble->setEnergy(); // set initial energy
    potE = ensemble->getPotEnergy();
    eFixed = ensemble->gettotEfixed();
    kNew = eFixed - (potE + lrc); // determine inital kinetic energy

    if (kNew > 0)
        ensemble->setKineticE(kNew);
    else {
        cout << "Initial configutaion has negative kinetic energy. Aborting." << endl;
        exit(0);
    }
    cout << "Results will be directed to the file \"mc_nve.out\"" << endl;

    for (i = 0; i < nStep; i++) {

        for (k = 0; k < num; ++k) {
            ++counter;
            potE = ensemble->getPotEnergy();
            kOld = ensemble->getKineticE();

            // Select atom randomly
            index = (int)num * random(&seed);
            position = atom[index]->getPosition();
            atom[index]->setTrialPosition(position);

            // determine energy at current position
            eOld = ensemble->getTrialPotE(index, atom);

            // find trial position
            trial[0] = position[0] + rMax * (2.0 * random(&seed) - 1.0);
            trial[1] = position[1] + rMax * (2.0 * random(&seed) - 1.0);
            trial[2] = position[2] + rMax * (2.0 * random(&seed) - 1.0);

            // apply periodic boundary conditions
            trial[0] -= length * nearestInt(trial[0], length);
            trial[1] -= length * nearestInt(trial[1], length);
            trial[2] -= length * nearestInt(trial[2], length);
            atom[index]->setTrialPosition(trial);

            // calculate energy at trial position
            eNew = ensemble->getTrialPotE(index, atom);
            eDelta = eNew - eOld;

            trialPotE = potE + eDelta;
            kNew = eFixed - (trialPotE + lrc);

            if (kNew > 0) {
                wFactor = pow(kNew / kOld, f);
                if (wFactor > 1 || wFactor > random(&seed)) // accept move
                {
                    atom[index]->reSetPosition();
                    ensemble->setKineticE(kNew);
                    ensemble->reSetEnergy(trialPotE);
                    accept++;
                }
            }
        }

        tally = accept / (double)counter;
        if (i < nEquil) {
            if (tally < 0.5)
                rMax *= 0.95;
            else
                rMax *= 1.05;
        }

        if (i >= nEquil) {
            potE = ensemble->getPotEnergy() + lrc;
            kinE = ensemble->getKineticE();
            totE = potE + kinE;
            temp = 2.0 * kinE / (3.0 * num);
            avPotEnergy += potE;
            avKinEnergy += kinE;
            avTotEnergy += totE;
            avTemperature += temp;

            if ((i != nEquil) && ((i % nSize) == 0)) {
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

    for (i = 0; i < n; i++) {
        erPotEnergy += pow((accumPotE[i] - avPotEnergy), 2);
        erKinEnergy += pow((accumKinE[i] - avKinEnergy), 2);
        erTotEnergy += pow((accumTotE[i] - avTotEnergy), 2);
        erTemperature += pow((accumTemp[i] - avTemperature), 2);
    }

    out << "Average Potential Energy/N:\t" << avPotEnergy / num << " +/- "
        << sqrt(erPotEnergy) / (num * nTotal) << endl;
    out << "Average Kinetic Energy/N:\t " << avKinEnergy / num << " +/- "
        << sqrt(erKinEnergy) / (num * nTotal) << endl;
    out << "Average Total Energy/N:\t\t" << avTotEnergy / num << "  +/- "
        << sqrt(erTotEnergy) / (num * nTotal) << endl;
    out << "Average Temperature:\t\t " << avTemperature << "  +/- " << sqrt(erTemperature) / nTotal
        << endl;
    out << "Acceptance Rate:\t\t " << 100 * tally << " %" << endl;
    out.close();
}

void MonteCarlo::readInNVE(int interPotential) {
    int i, j;
    double *molFract, **epsilon, **sigma, **rCut, potE, density;
    int dimensions = 3, numComp, numAtom, *comp;
    Atom *newAtom, **atoms;

    std::string configPath = std::string(DATA_PATH) + "/config_mc.toml";

    try {
        auto config = toml::parse_file(configPath);

        // Read number of components
        numComp = config["ensemble"]["components"].value_or(0);
        if (numComp <= 0) {
            std::cerr << "Number of components must be > 0!" << std::endl;
            return;
        }

        // Allocate memory for components
        comp = new int[numComp];

        // Read number of atoms
        numAtom = config["ensemble"]["atoms"].value_or(0);
        if (numAtom < 4) {
            std::cerr << "At least 4 atoms are required!" << std::endl;
            return;
        }

        // Allocate memory and read mole_fraction
        molFract = new double[numComp];
        double moleFractionValue = config["ensemble"]["mole_fraction"].value_or(0.0);
        if (moleFractionValue <= 0.0 || moleFractionValue > 1.0) {
            std::cerr << "Invalid mole fraction value!" << std::endl;
            return;
        }
        for (i = 0; i < numComp; i++) {
            molFract[i] = moleFractionValue; // Use the same value for all components
        }

        // Allocate memory and read mass
        double* mass = new double[numComp];
        double massValue = config["ensemble"]["mass"].value_or(0.0);
        if (massValue <= 0.0) {
            std::cerr << "Invalid mass value!" << std::endl;
            return;
        }
        for (i = 0; i < numComp; i++) {
            mass[i] = massValue; // Use the same value for all components
        }

        // Read total energy and density
        potE = config["results"]["total_energy"].value_or(0.0);
        density = config["results"]["density"].value_or(0.0);
        potE *= numAtom;

        if (interPotential == 1) { // Lennard-Jones potential
            double epsilonValue = config["results"]["epsilon"].value_or(0.0);
            double sigmaValue = config["results"]["sigma"].value_or(0.0);
            double rCutValue = config["results"]["rcut"].value_or(0.0);

            if (epsilonValue <= 0.0 || sigmaValue <= 0.0 || rCutValue <= 0.0) {
                std::cerr << "Invalid Lennard-Jones parameters!" << std::endl;
                return;
            }

            // Allocate memory for interaction parameters
            epsilon = getMemory(numComp);
            sigma = getMemory(numComp);
            rCut = getMemory(numComp);

            for (i = 0; i < numComp; i++) {
                for (j = 0; j < numComp; j++) {
                    epsilon[i][j] = epsilonValue;
                    sigma[i][j] = sigmaValue;
                    rCut[i][j] = rCutValue;

                    // Adjust for indistinguishable pairs
                    if (j != i) {
                        epsilon[j][i] = epsilon[i][j];
                        sigma[j][i] = sigma[i][j];
                        rCut[j][i] = rCut[i][j];
                    }
                }
            }

            // Allocate and initialize atoms
            atoms = new Atom*[numAtom];
            int* atnum = new int[numComp + 1];
            atnum[0] = 0;

            // Calculate the number of atoms of each type from the mole fraction
            for (i = 0; i < numComp; i++) {
                double numType = numAtom * molFract[i];
                int rem = static_cast<int>(numType);
                atnum[i + 1] = (numType - rem >= 0.5) ? rem + 1 : rem;
            }

            while (atnum[numComp] < numAtom) atnum[numComp]++;
            for (i = 0; i < numComp; i++) atnum[i + 1] += atnum[i];

            for (i = 0; i < numComp; i++) {
                for (int j = atnum[i]; j < atnum[i + 1]; j++) {
                    newAtom = new LJatom(i, mass[i], epsilon, sigma, rCut, dimensions);
                    atoms[j] = newAtom;
                }
            }
        }

        // Create the NVE ensemble
        ensemble = new NVEensemble(atoms, numComp, numAtom, potE, density, molFract, comp);

    } catch (const toml::parse_error& err) {
        std::cerr << "Error parsing TOML file: " << err << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Unexpected error: " << e.what() << std::endl;
    }
}

// Method: getnSize()
// Usage: n = getSize();
// --------------------
// Return the size of the blocks for averaging
// ensemble properties.

int MonteCarlo::getnSize() { return nSize; }

// Method: getnEquil()
// Usage: n = getnEquil();
// -----------------------
// Return the number of steps prior to
// equilibrium.

int MonteCarlo::getnEquil() { return nEquil; }

// Method: getnStep()
// Usage: n = getnStep();
// ----------------------
// Return the total number of simulation
// steps.

int MonteCarlo::getnStep() { return nStep; }

// Method:  run
// Usage:  run();
// --------------
// Runs the appropriate simularion run method
// determined by the choice of ensemble.

void MonteCarlo::run() {
    switch (theEnsemble) {
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
/*
MonteCarlo::MonteCarlo() {
    int potential;
    // ifstream in;
    // in.open("mc.dat");

    // Use the DATA_PATH macro defined in CMake
    string dataPath = string(DATA_PATH) + "/mc.dat";

    cout << "Attempting to open: " << dataPath << endl;

    ifstream in(dataPath);
    if (in.fail()) {
        cout << "Cannot open " << dataPath << "!\n";
        return;
    }

    if (in.fail()) {
        cout << "Unable to open mc.dat for MC parameters" << endl;
        return;
    }

    in >> nStep;
    in >> nEquil;
    in >> nSize;
    in >> theEnsemble >> potential;
    in.close();

    switch (theEnsemble) {
    case 1: // NVE ensemble
        readInNVE(potential);
        break;
        // other ensembles here
    }
}
*/

MonteCarlo::MonteCarlo() {
    int potential;
    // Path to the TOML configuration file
    std::string configPath = std::string(DATA_PATH) + "/config_mc.toml";

    std::cout << "Attempting to parse: " << configPath << std::endl;

    try {
        // Parse the TOML file
        auto config = toml::parse_file(configPath);

        // Read Monte Carlo parameters
        nStep = config["monte_carlo"]["total_cycles"].value_or(0);
        nEquil = config["monte_carlo"]["equilibration_cycles"].value_or(0);
        nSize = config["monte_carlo"]["block_size"].value_or(0);

        // Read ensemble and potential type
        std::string ensembleType = config["ensemble"]["type"].value_or("");
        std::string potentialType = config["ensemble"]["potential"].value_or("");

        if (ensembleType.empty() || potentialType.empty()) {
            std::cerr << "Missing or invalid ensemble or potential type in configuration!" << std::endl;
            return;
        }

        if (ensembleType == "NVE" && potentialType == "Lennard-Jones") {
            theEnsemble = 1; // Map "NVE" to 1  (global... FIXME)
            potential = 1;   // Map "Lennard-Jones" to 1
            readInNVE(potential); // Read NVE-specific settings
        } else {
            std::cerr << "Unsupported ensemble or potential type: "
                      << "Ensemble = " << ensembleType
                      << ", Potential = " << potentialType << std::endl;
        }
    } catch (const toml::parse_error& err) {
        // Handle TOML parsing errors
        std::cerr << "Error parsing TOML file: " << err << std::endl;
    } catch (const std::exception& e) {
        // Handle other unexpected errors
        std::cerr << "Unexpected error: " << e.what() << std::endl;
    }
}
