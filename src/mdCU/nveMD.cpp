// File: NVEensembleMD.cpp
// -------------------------
// File containing functions to implement the
// NVEensemble class.

// This code was specifically developed to illustrate concepts in the accompanying book:
// R. J. Sadus, "Molecular Simulation of Fluids: Theory, Algorithms, Object-Orientation,
// and Parallel Computing," 2nd Ed. (Elsevier, Amsterdam, 2023). It can be used freely for
// any not for profit purpose or academic research application. The code has been validated,
// but it would be nonetheless prudent to test it further before publishing any results.
// Check the book's website for any subsequent updates.

#include "nveMD.h"
#include "auxfunc.h"
#include "ljforceMD.h"
#include <iostream>
#include <math.h>

using namespace std;

// Method: setKineticE
// Usage: setKineticE();
// ----------------------
// Determines the kinetic energy from the velocities.

void NVEensemble::setKineticE() {
    int i;
    double sum = 0.0, *velocity;

    for (i = 0; i < numAtom; i++) {
        velocity = atoms[i]->getVelocity();
        sum += atoms[i]->getMass() *
               (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2]);
    }

    kineticE = sum / 2;
}

// Method: setKineticE
// Usage:  setKineticE(kinE)
// _________________________
// Sets kinetic energy after Gear Corrector step

void NVEensemble::setKineticE(double kinE) { kineticE = kinE; }

// Method: scaleVelocity
// Usage: scaleVelocity();
// ---------------------------
// Sacles the velocities during the equilibration phase
// of the simulation. The scaling ensures that the simulation
// is performed at the specified temperature. The kinetic
// energy must be determined at least once before it is invoked.

void NVEensemble::scaleVelocity() {
    int i;
    double tScale, mass, kinE, temp, sumVSq = 0.0, *velocity;

    kinE = getKineticE();
    temp = 2 * kinE / (3 * numAtom - 3);

    tScale = sqrt(getTemp() / temp);

    // scale velocities to reproduce the desired
    // temperature (temp) and acumulate summations
    // for momentum adjustment

    for (i = 0; i < numAtom; i++) {
        velocity = atoms[i]->getVelocity();

        velocity[0] *= tScale;
        velocity[1] *= tScale;
        velocity[2] *= tScale;

        atoms[i]->setVelocity(velocity);
        ;
    }
}

// Method: initialVelocity
// Usage: initialVelocity();
// -------------------------
// Invokes a random number generator to assign initial
// velocities to molecules at the start of a molecular dynamics
// simulation.

void NVEensemble::initialVelocity() {
    int seed; // seed for random number generator
    int i;
    double vX, vY, vZ, vSq, rValue;
    double sumX = 0.0, sumY = 0.0;
    double sumZ = 0.0;
    double* velocity;

    // assign random velocities to molecules in the range of
    // 0 to 1

    seed = -29;

    for (i = 0; i < numAtom; i++) {
        vX = random(&seed);
        vY = random(&seed);
        vZ = random(&seed);

        vSq = vX * vX + vY * vY + vZ * vZ;

        rValue = sqrt(vSq);

        vX /= rValue;
        vY /= rValue;
        vZ /= rValue;

        sumX += vX;
        sumY += vY;
        sumZ += vZ;

        velocity = new double[3]; // create new array with 3 elements
        velocity[0] = vX;         // assign x-component of velocity
        velocity[1] = vY;         // assign y-component of velocity
        velocity[2] = vZ;         // assign z-component of velocity

        atoms[i]->setVelocity(velocity); // adjust velocity of atom i
    }
    // adjust velocities so that the total linear momentum
    // is zero

    for (i = 0; i < numAtom; i++) {
        velocity = atoms[i]->getVelocity();
        velocity[0] -= sumX / ((double)numAtom);
        velocity[1] -= sumY / ((double)numAtom);
        velocity[2] -= sumZ / ((double)numAtom);
        atoms[i]->setVelocity(velocity);
    }
}

// Method: getEnergyLRC
// Usage: n = getEnergyLRC();
// --------------------------
// Gets the long range energy correction.

double NVEensemble::getEnergyLRC() { return energyLRC; }

// Method: getVirialLRC
// Usage: n = getVirialLRC();
// --------------------------
// Gets the long range virial correction.

double NVEensemble::getVirialLRC() { return virialLRC; }

// Method: initialCoord
// Usage: initialCoord()
// ---------------------
// Places atoms upon a lattice to get initial coordinates

void NVEensemble::initialCoord() {
    int i, j, x, y, z, offset;
    double cells, dcells, cellL, halfCellL, *tempPos, *position;

    // Determine the number of unit cells in each coordinate
    // direction

    dcells = pow(0.25 * (double)numAtom, 1.0 / 3.0);
    cells = (int)nearestInt(dcells, 1.0);

    // check if numAtom is an non-fcc number of molecules
    // and increase the number of cells if necessary

    while ((4 * cells * cells * cells) < numAtom)
        cells = cells + 1;

    // Determine length of the unit cell

    cellL = boxLen / (double)cells;
    halfCellL = 0.5 * cellL;

    // Construct the unit cell
    // point to atoms position
    position = atoms[0]->getPosition();
    position[0] = 0.0;
    position[1] = 0.0;
    position[2] = 0.0;

    position = atoms[1]->getPosition();
    position[0] = halfCellL;
    position[1] = halfCellL;
    position[2] = 0.0;

    position = atoms[2]->getPosition();
    position[0] = 0.0;
    position[1] = halfCellL;
    position[2] = halfCellL;

    position = atoms[3]->getPosition();
    position[0] = halfCellL;
    position[1] = 0.0;
    position[2] = halfCellL;

    for (i = 4; i < numAtom; i++) {
        position = atoms[i]->getPosition();
        position[0] = 0.0;
        position[1] = 0.0;
        position[2] = 0.0;
    } // init all other atoms to 0

    // Build the lattice from the unit cell by
    // repeatly translating the four vectors of
    // the unit cell through a distance cellL in
    // the x, y and z directions

    offset = 0;

    for (z = 1; z <= cells; z++)
        for (y = 1; y <= cells; y++)
            for (x = 1; x <= cells; x++) {
                for (i = 0; i < 4; i++) {
                    j = i + offset;
                    if (j < numAtom) {
                        tempPos = atoms[j]->getPosition();
                        position = atoms[i]->getPosition();
                        tempPos[0] = position[0] + cellL * (x - 1);
                        tempPos[1] = position[1] + cellL * (y - 1);
                        tempPos[2] = position[2] + cellL * (z - 1);
                    }
                }

                offset = offset + 4;
            }

    // Shift centre of box to the origin.

    for (i = 0; i < numAtom; i++) {
        tempPos = atoms[i]->getPosition();
        tempPos[0] -= halfCellL;
        tempPos[1] -= halfCellL;
        tempPos[2] -= halfCellL;
        atoms[i]->setPosition(tempPos);
        ;
    }
}

// Method: setforce
// Usage: setForce();
// ensemble.cpp
// ------------------
// initiate the force calculations using LJforce

void NVEensemble::setForce() { theForce->setForce(numAtom, boxLen, &potEnergy, &virial); }

// Method: lrc
// Usage: lrc();
// -------------
// initiate the force calculations for long range corrections

void NVEensemble::lrc() { theForce->lrc(numComp, comp, boxVol, &energyLRC, &virialLRC); }

// Method: NVEensemble
// Usage: NVEensemble;
// -------------------
// Instantiate the NVEensemble class

NVEensemble::NVEensemble(Atom** atoms, int numComp, int numAtom, double temperature, double density,
                         double* molFract, int* comp)
    : Ensemble(atoms, numComp, numAtom, temperature, density, molFract, comp) {
    setVolume();
    setLength();
    setComp();
    initialCoord();
    theForce = new LJforce(atoms);
}
