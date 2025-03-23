// File: nveMC.cpp
// ----------------
// File containing functions to implement the
// NVEensemble class.

// This code was specifically developed to illustrate concepts in the accompanying book:
// R. J. Sadus, "Molecular Simulation of Fluids: Theory, Algorithms, Object-Orientation,
// and Parallel Computing," 2nd Ed. (Elsevier, Amsterdam, 2023). It can be used freely for
// any not for profit purpose or academic research application. The code has been validated,
// but it would be nonetheless prudent to test it further before publishing any results.
// Check the book's website for any subsequent updates.

#include "nveMC.h"
#include "auxfunc.h"
#include "ljenergyMC.h"
#include <math.h>
#include <string.h>

// Method: getEnergyLRC
// Usage: n = getEnergyLRC();
// --------------------------
// Gets the long range energy correction.

double NVEensemble::getEnergyLRC() { return energyLRC; }

// Method: initialCoord
// Usage: initialCoord()
// ---------------------
// Places atoms upon a lattice to get initial coordinates

void NVEensemble::initialCoord() {
    int i, j, x, y, z, offset;
    double cells, dcells, cellL, halfCellL, *tempPos, *position;
    double(*holdPos)[3] = new double[numAtom][3];
    tempPos = new double[numAtom];
    position = new double[numAtom];

    // Determine the number of unit cells in each coordinate
    // direction

    dcells = pow(0.25 * numAtom, 1.0 / 3);
    cells = (int)nearestInt(dcells, 1.0);

    // Check if numAtom is an non-fcc number of molecules
    // and increase the number of cells if necessary

    while ((4 * cells * cells * cells) < numAtom)
        cells = cells + 1;

    // Determine length of the unit cell

    cellL = boxLen / (double)cells;
    halfCellL = cellL / 2;

    // Construct the unit cell
    // Point to atoms position
    position = atoms[0]->getPosition();
    position[0] = 0.0;
    position[1] = 0.0;
    position[2] = 0.0;
    atoms[0]->setPosition(position);

    position = atoms[1]->getPosition();
    position[0] = halfCellL;
    position[1] = halfCellL;
    position[2] = 0.0;
    atoms[1]->setPosition(position);

    position = atoms[2]->getPosition();
    position[0] = 0.0;
    position[1] = halfCellL;
    position[2] = halfCellL;
    atoms[2]->setPosition(position);

    position = atoms[3]->getPosition();
    position[0] = halfCellL;
    position[1] = 0.0;
    position[2] = halfCellL;
    atoms[3]->setPosition(position);

    for (i = 4; i < numAtom; i++) {
        position = atoms[i]->getPosition();
        position[0] = 0.0;
        position[1] = 0.0;
        position[2] = 0.0;
        atoms[i]->setPosition(position);

    } // initalise all other atoms to 0

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
                        position = atoms[i]->getPosition();
                        holdPos[j][0] = position[0] + cellL * (x - 1);
                        holdPos[j][1] = position[1] + cellL * (y - 1);
                        holdPos[j][2] = position[2] + cellL * (z - 1);
                    }
                }

                offset = offset + 4;
            }

    // Shift centre of box to the origin.

    for (i = 0; i < numAtom; i++) {
        tempPos = atoms[i]->getPosition();
        tempPos[0] = holdPos[i][0] - halfCellL;
        tempPos[1] = holdPos[i][1] - halfCellL;
        tempPos[2] = holdPos[i][2] - halfCellL;
        atoms[i]->setPosition(tempPos);
    }
}
// Method: reSetEnergy
// Usage:  reSetEnergy(energy);
// ----------------------------
// Re-sets the energy of the ensemble following
// an acceted move

void NVEensemble::reSetEnergy(double eng) { potEnergy = eng; }

// Method: setEnergy
// Usage: setEnergy();
// -------------------
// Set the energy calculations using LJenergy.
// This is an N*N calculation to set the energy
// at the start of the simulation.

void NVEensemble::setEnergy() { theEnergy->setEnergy(numAtom, boxLen, &potEnergy); }

// Method: getTrialPotE
// Usage: n = getTrialPotE(index, atoms)
// ---------------------------------
// Performs N calculations to update energy before and after
// index atom is moved.

double NVEensemble::getTrialPotE(int index, Atom** a) {
    return theEnergy->getTrialPotE(numAtom, boxLen, index, a);
}
// Method: lrc
// Usage: lrc();
//  ------------
// initiate the force calculations for long range corrections

void NVEensemble::lrc() { theEnergy->lrc(numComp, comp, boxVol, &energyLRC); }

// Method: NVEensemble
// Usage: NVEensemble;
// -------------------
//  instantiate the NVEensemble class

NVEensemble::NVEensemble(Atom** atoms, int numComp, int numAtom, double totEfixed, double density,
                         double* molFract, int* comp)
    : Ensemble(atoms, numComp, numAtom, totEfixed, density, molFract, comp) {
    setVolume();
    setLength();
    setComp();
    initialCoord();
    theEnergy = new LJenergy(atoms);
}
