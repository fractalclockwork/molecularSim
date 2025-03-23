// File: energyMC.cpp
// ----------------
// File containing functions to implement the abstract
// Energy class.

// This code was specifically developed to illustrate concepts in the accompanying book:
// R. J. Sadus, "Molecular Simulation of Fluids: Theory, Algorithms, Object-Orientation,
// and Parallel Computing," 2nd Ed. (Elsevier, Amsterdam, 2023). It can be used freely for
// any not for profit purpose or academic research application. The code has been validated,
// but it would be nonetheless prudent to test it further before publishing any results.
// Check the book's website for any subsequent updates.

#include "energyMC.h"

// constructor
// Method: Energy
// Usage:   n = Energy(atoms);
// --------------------------
// Used to construct the force component of any derived energy classes.

Energy::Energy(Atom** theAtoms) { atoms = theAtoms; }

// destructor
Energy::~Energy() {}
