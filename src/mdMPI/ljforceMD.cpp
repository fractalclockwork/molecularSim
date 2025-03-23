// File: LJforceMD.cpp
// -------------------
// File containing methods to implement the LJforce
// class.

// This code was specifically developed to illustrate concepts in the accompanying book:
// R. J. Sadus, "Molecular Simulation of Fluids: Theory, Algorithms, Object-Orientation,
// and Parallel Computing," 2nd Ed. (Elsevier, Amsterdam, 2023). It can be used freely for
// any not for profit purpose or academic research application. The code has been validated,
// but it would be nonetheless prudent to test it further before publishing any results.
// Check the book's website for any subsequent updates.

#include "ljforceMD.h"
#include "auxfunc.h"
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <string.h>

using namespace std;

// Method: setForce
// Usage: setForce();
// ------------------
// The ljForce function calculates the force experieced
// by all num atoms due to pairwise interatomic interaction.
// The forces are calculated using the Lennard-Jones (6-12)
// potential. The potential is truncated at a distance
// rCut and long range corrections must be applied outside
// the function to obtain the full contributions to both
// energy and the virial.

void LJforce::setForce(int num, double length, double* potEnergy, double* virial) {
    int i, j, k, l, kindi, kindj;
    double rXi, rYi, rZi;           // position vectors for atom i
    double rXj, rYj, rZj;           // position vectors for atom j
    double rXij, rYij, rZij;        // pair seperation vectors
    double rijSq;                   // interatomic seperation squared
    double fXi, fYi, fZi;           // force vectors for atom i
    double fij, fXij, fYij, fZij;   // force between atoms i and j
    double rCutSq, sigmaSq;         // squared cut and sigma values
    double sigma2, sigma6, sigma12; // multiples of LJ sigma
    double pot;                     // Potential between 2 atoms
    double* tempPos;                // postion vectors
    double* tempPosj;               // postion vector
    double potECal, virialCal;

    // Additonal variables for MPI implementation
    int myNode, totalNodes, lastNode;
    int start, end, numCal;
    double potECalR, virialCalR; // Receiving variables for MPI_Allreduce
    double* f = new double[3 * num];

    // Zero the force  arrays.
    for (i = 0; i < num; i++) {
        k = 3 * i;
        f[k] = f[k + 1] = f[k + 2] = 0;
    }

    potECal = 0.0;
    virialCal = 0.0;

    // Set-up MPI
    MPI_Comm_size(MPI_COMM_WORLD, &totalNodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myNode);

    lastNode = totalNodes - 1; // Node count starts at 0
    numCal = num / totalNodes; // Distribute calculation evenly between nodes
    start = myNode * numCal;   // Different starting point for each node

    if (myNode == lastNode) // Different end point for each node
        end = num;
    else
        end = (myNode + 1) * numCal;

    int step = 1;

    // An alternative MPI implementation can be obtained by setting:
    // step = totalNodes;
    // start = myNode;
    // end = num;

    // Calculate forces experienced by atoms due
    // interatomic pair interactions.

    for (i = start; i < end; i += step) {
        k = 3 * i;
        // Select molecule i
        kindi = atoms[i]->getType();

        rXi = atoms[0]->r[k];
        rYi = atoms[0]->r[k + 1];
        rZi = atoms[0]->r[k + 2];

        fXi = f[k];
        fYi = f[k + 1];
        fZi = f[k + 2];

        for (j = 0; j < num; j++) {
            if (j != i) {
                kindj = atoms[j]->getType();

                l = 3 * j;
                rXj = atoms[0]->r[l];
                rYj = atoms[0]->r[l + 1];
                rZj = atoms[0]->r[l + 2];

                // Calculate pair seperation
                rXij = rXi - rXj;
                rYij = rYi - rYj;
                rZij = rZi - rZj;

                // Apply periodic boundary
                rXij -= length * nearestInt(rXij, length);
                rYij -= length * nearestInt(rYij, length);
                rZij -= length * nearestInt(rZij, length);

                rijSq = rXij * rXij + rYij * rYij + rZij * rZij;

                // Calculate forces for seprations
                // below the cutoff
                rCutSq = rCut[kindi][kindj] * rCut[kindi][kindj];

                if (rijSq <= rCutSq) {
                    // Determine potential between i and j

                    sigmaSq = sigma[kindi][kindj] * sigma[kindi][kindj];
                    sigma2 = sigmaSq / rijSq;
                    sigma6 = sigma2 * sigma2 * sigma2;
                    sigma12 = sigma6 * sigma6;
                    pot = sigma12 - sigma6;

                    potECal += 4 * epsilon[kindi][kindj] * pot;
                    virialCal += 24 * epsilon[kindi][kindj] * (pot + sigma12) / 3;

                    // Calculate forces between atoms
                    fij = 24 * epsilon[kindi][kindj] * (pot + sigma12) / rijSq;
                    fXij = fij * rXij;
                    fYij = fij * rYij;
                    fZij = fij * rZij;
                    fXi += fXij;
                    fYi += fYij;
                    fZi += fZij;
                }
            }
        }

        f[k] = fXi;
        f[k + 1] = fYi;
        f[k + 2] = fZi;
    }

    // MPI Reductions

    MPI_Allreduce(&potECal, &potECalR, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&virialCal, &virialCalR, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    MPI_Allreduce(f, atoms[0]->f, 3 * num, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    delete[] f;

    *potEnergy = potECalR / 2;
    *virial = virialCalR / 2;
}

// Method: lrc
// Usage: n = lrc(num, nComp, volume, energyLRC, virialLRC);
// ---------------------------------------------------------
// ljLRC calculates the long range correction terms to both
// the energy and virial for an ensemble interacting via the
// lennard-Jones (12-6) potential.

void LJforce::lrc(int num, int* nComp, double boxV, double* energyLRC, double* virialLRC) {
    int i, j;
    const double PI = 3.141592654;
    double sigmaSq, sig3, rCutSq, rCut3, sig3R, sig9R, den;

    *energyLRC = 0.0;
    *virialLRC = 0.0;

    for (i = 0; i < num; i++)
        for (j = 0; j < num; j++) {
            sigmaSq = sigma[i][j] * sigma[i][j];
            rCutSq = rCut[i][j] * rCut[i][j];
            sig3 = sigma[i][j] * sigmaSq;
            rCut3 = rCut[i][j] * rCutSq;
            sig3R = sig3 / rCut3;
            sig9R = sig3R * sig3R * sig3R;
            den = nComp[i] * nComp[j] * sig3 / boxV;
            *energyLRC += (8 / 9.0) * den * PI * epsilon[i][j] * (sig9R - 3 * sig3R);
            *virialLRC += (16 / 9.0) * den * PI * epsilon[i][j] * (2 * sig9R - 3 * sig3R);
        }
}

// constructor
LJforce::LJforce(Atom** theAtoms) : Force(theAtoms) {
    epsilon = atoms[0]->getEpsilon();
    sigma = atoms[0]->getSigma();
    rCut = atoms[0]->getrCutOff();
}
