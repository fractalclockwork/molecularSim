// File: mainMD.cpp
// _________________
// File containing the main() function for the MD program.

// This code was specifically developed to illustrate concepts in the accompanying book:
// R. J. Sadus, "Molecular Simulation of Fluids: Theory, Algorithms, Object-Orientation,
// and Parallel Computing," 2nd Ed. (Elsevier, Amsterdam, 2023). It can be used freely for
// any not for profit purpose or academic research application. The code has been validated,
// but it would be nonetheless prudent to test it further before publishing any results. 
// Check the book's website for any subsequent updates.

#include "simMD.h"

int main(int argc, char **argv)
{
    
 Simulation sim;
 sim.readSimParameters();

 return 0;
}
