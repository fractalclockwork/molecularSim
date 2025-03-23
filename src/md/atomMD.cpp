// File: atomMD.cpp
// ---------------- 
// File containing the method bodies to implement the
// atom class.
 
// This code was specifically developed to illustrate concepts in the accompanying book:
// R. J. Sadus, "Molecular Simulation of Fluids: Theory, Algorithms, Object-Orientation,
// and Parallel Computing," 2nd Ed. (Elsevier, Amsterdam, 2023). It can be used freely for
// any not for profit purpose or academic research application. The code has been validated,
// but it would be nonetheless prudent to test it further before publishing any results. 
// Check the book's website for any subsequent updates.

#include "atomMD.h"
#include "auxfunc.h"
#include <iostream>
#include <fstream>
#include <math.h>

// Method:  setType
// Usage:   setType(int atomType);
// -------------------------------
//  Used to identify the atom as belonging to a certain
//  _type_ of component within the ensemble

void Atom::setType(int t)
{
  type = t;
}

// Method: setrCutOff
// Usage:  setrCutOff(rCutOffValue);
// -------------------------------------------------
// Used to assign the values for the rCutOffs for the various
// atom types that make up the ensemble

void Atom::setrCutOff(double **newrCut)
{
    rCutOff = newrCut;
}

// Method:  setMass
// Usage:   setMass(double atomicMass);
// ------------------------------------
// Assign the mass of the atom

void Atom::setMass(double newMass)
{
  mass = newMass;
}

// Method:  setAcceleration
// Usage:   setAcceleration(atomicAcceleration);
// -----------------------------------------------------
// Sets the pointer within the atom aboject to point to a
// new 1D array of doubles, representing the acceleration
// of the atom in its various dimensions.  Assuming a 3D
// simulation, the required array would have (only) three
// elements, where:
// atomicAcceleration[0] = acceleration in X dimension
// atomicAcceleration[1] = acceleration in Y dimension
// atomicAcceleration[2] = acceleration in Z dimension
// Note:  this method requires an _array_ to be passed to it
// which has already had memory allocated to it. 

void Atom::setAcceleration(double *newAccel)
{
  for(int i = 0; i < 3; i++)
    acceleration[i] = newAccel[i];
}

// Method:  setPosition
// Usage:   setPosition(atomicPosition);
// ---------------------------------------------
// Sets the pointer within the atom aboject to point to a
// new 1D array of doubles, representing the position
// of the atom in its various dimensions.  Assuming a 3D
// simulation, the required array would have (only) three
// elements, where:
// atomicPosition[0] = position in X dimension
// atomicPosition[1] = position in Y dimension
// atomicPosition[2] = position in Z dimension
// Note:  this method requires an _array_ to be passed to it
// which has already had memory allocated to it.
	
void Atom::setPosition(double *newPos)
{
  for(int i = 0; i < 3; i++)
    position[i] = newPos[i];
}

// Method:  setVelocity
// Usage:   setVelocity(atomicVelocity);
// ---------------------------------------------
// Sets the pointer within the atom aboject to point to a
// new 1D array of doubles, representing the velocity
// of the atom in its various dimensions.  Assuming a 3D
// simulation, the required array would have (only) three
// elements, where:
// atomicVelocity[0] = velocity in X dimension
// atomicVelocity[1] = velocity in Y dimension
// atomicVelocity[2] = velocity in Z dimension
// Note:  this method requires an _array_ to be passed to it
// which has already had memory allocated to it. 

void Atom::setVelocity(double *newV)
{
 for(int i = 0; i < 3; i++)
    velocity[i] = newV[i];
}

// Method:  setHighTimeDerivs
// Usage:   setHighTimeDerivs(timeDerivatives);
// ---------------------------------------------
// Sets the pointer within the atom aboject to point to a
// new 2D array of doubles, representing the derivatives
// of time after acceleration and velocity fro that atom,
// in its various dimensions.  Assuming a 3D simulation,
// the required array would be 2D where:
// the first dimension represents the derivative of time, and
// the second dmiension represents the dimension in space.
// For example,
// atomicTimeStep[x][y] refers to the derivative of time for:
// - derivative x, where x is 0 for the third derivative of time,
// 1 is the fourth derivative of time, and 2 is the fifth, and
// - dimension y, where y is 0 for the X dimension, 1 for the Y
// dimension, and 2 for the Z dimension
// Note:  this method requires an _array_ to be passed to it
// which has already had memory allocated to it.  

void Atom::setHighTimeDerivs(double **newHighTimeDerivs)
{
 for(int i = 0; i < 3; i++)
   for(int j = 0; j < 3; j++) 
     highTimeDerivs[i][j] = newHighTimeDerivs[i][j];
}

// Method:  setForce
// Usage:   setForce(double *atomicForce);
// ---------------------------------------------
// Sets the pointer within the atom object to point to a
// new 1D array of doubles, representing the force acting
// on the atom in its various dimensions.  Assuming a 3D
// simulation, the required array would have (only) three
// elements, where:
// atomicForce[0] = force in X dimension
// atomicForce[1] = force in Y dimension
// atomicForce[2] = force in Z dimension
// Note:  this method requires an _array_ to be passed to it
// which has already had memory allocated to it.

void Atom::setForce(double *newForce)
{
    for(int i = 0; i < 3; i++)
     force[i] = newForce[i];
}

void Atom::setForce(double fX, double fY, double fZ)
{
  force[0] = fX;
  force[1] = fY;
  force[2] = fZ;
}

void Atom::setForce(int i, double *newForce)
{
   for(int i = 0; i < 3; i++)
      forceOld[i] = newForce[i];
}

// Access (Get) Methods
// Usage:   n = getXxxx();
// --------------------
// These methods return the _value_ stored for the mass or
// type of the atom.

// Method: getType
// Usage: p = getType();
// ---------------------- 
// Gets the type associated with each atom.

int Atom::getType()
{
  return type;
}

// Method: getMass
// Usage: n = getMass();
// --------------------- 
// Get values of the atomic masses.

double Atom::getMass()
{
  return mass;
}

// Access (Get) Methods
// Usage:   n = getXxxxx();
// ------------------------
// These methods return a reference to the arrays which hold
// the values for rCutOff, acceleration, position, velocity,
// time step, and force, as described above their respective
// set method.

// Method: getrCutOff
// Usage: n = getrCutOff();
// ---------------------
// Returns the cut off distances for atom pairs

double **Atom::getrCutOff()
{
  return rCutOff;
}

// Method: getAcceleration
// Usage: n = getAcceleration();
// ---------------------------- 
// Return acceleration of atom

double *Atom::getAcceleration()
{
  return acceleration;
}

// Method: getPosition
// Usage: n = getPosition();
// -------------------------
// Return position of atom

double *Atom::getPosition()
{
  return position;
}

// Method: getVelocity
// Usage: n = getVelosity();
// ------------------------- 
// Return the velocity of the atom

double *Atom::getVelocity()
{
  return velocity;
}

// Method:  getHighTimeDerivs
// Usage: n = getHighTimeDerivs();
// ------------------------------- 
// Return the derivatives of time

double **Atom::getHighTimeDerivs()
{
  return highTimeDerivs;
}

// Method:  getForce
// Usage: n = getForce();
// ----------------------
// Returns reference to the force array in the atom object

double *Atom:: getForce()
{
  return force;
}

double *Atom:: getForce(int i)
{
  return forceOld;
}


// Constructor:Atom
// Usage: Atom atom;
// ----------------- 
// Builds the Atom class

Atom::Atom(int theType, double theMass, int dimensions, int derivatives)
{
  int i;
  mass = theMass;
  type = theType;
	
  //allocate memory for the arrays
  position = new double[dimensions];
  velocity = new double[dimensions];
  acceleration = new double[dimensions];
  force = new double[dimensions];
  forceOld = new double[dimensions];;
  highTimeDerivs = getMemory(derivatives, dimensions);
 
  //initialise values for acceleration and highTimeDerivs
  for (i = 0; i < 3; i++)
  {
    acceleration[i] = 0.0;
  }

  for(i = 0; i < 3; i++)
  {
     for (int j = 0; j < 3; j++)
     {
       highTimeDerivs[i][j] = 0.0;
      }
   }
}

// 2nd constructor for array reference construction
Atom::Atom(){}

