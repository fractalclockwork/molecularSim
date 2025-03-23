// File: gpcMD.cpp
// Class: GearPC

// This code was specifically developed to illustrate concepts in the accompanying book:
// R. J. Sadus, "Molecular Simulation of Fluids: Theory, Algorithms, Object-Orientation,
// and Parallel Computing," 2nd Ed. (Elsevier, Amsterdam, 2023). It can be used freely for
// any not for profit purpose or academic research application. The code has been validated,
// but it would be nonetheless prudent to test it further before publishing any results. 
// Check the book's website for any subsequent updates.

#include "gpcMD.h"
#include "gpcConstantsMD.h"
#include "auxfunc.h"
#include <iostream>

using namespace std;

// Method: gearPredict
// Usage: gearPredict(num, length);
// -------------------------------- 
// The function gearPredict uses Gear's 5th-order algorithm
// to predict new coordinates, accelerations, velocities,
// third and fourth time derivatives of position.

void GearPC::gearPredict(int num, double length, double delT)
{
  int i, j;
  double *position, *velocity, *acceleration, **highTimeDerivs;
  double c1, c2, c3, c4, c5;

  c1 = delT;
  c2 = c1*delT/2;
  c3 = c2*delT/3;
  c4 = c3*delT/4;
  c5 = c4*delT/5;
 
  for(i = 0; i < num; i++)
  {
    velocity = atoms[i]->getVelocity();
    acceleration = atoms[i]->getAcceleration();
    highTimeDerivs = atoms[i]->getHighTimeDerivs();

    for(j = 0; j < 3; j++)
    {
     //compoents in the x, y and z directions
     atoms[0]->r[3*i + j] += c1*velocity[j] + c2*acceleration[j]
		+  c3*highTimeDerivs[0][j]   //3rd derivative
		+  c4*highTimeDerivs[1][j]   //4th derivative
		+  c5*highTimeDerivs[2][j];  //5th derivative
     atoms[0]->r[3*i + j] -= length * nearestInt(atoms[0]->r[3*i + j], length);
     velocity[j] += c1*acceleration[j]
	        +  c2*highTimeDerivs[0][j] //3rd derivative
	        +  c3*highTimeDerivs[1][j] //4th derivative
	        +  c4*highTimeDerivs[2][j];//5th derivative
     acceleration[j] += c1*highTimeDerivs[0][j] //3rd derivative
		    +  c2*highTimeDerivs[1][j] //4th derivative
		    +  c3*highTimeDerivs[2][j];//5th derivative
     highTimeDerivs[0][j] += c1*highTimeDerivs[1][j] //4th derivative
	                  +  c2* highTimeDerivs[2][j]; //5th derivative
     highTimeDerivs[1][j]  += c1*highTimeDerivs[2][j]; //5th derivative
   }
   
   atoms[i]->setAcceleration(acceleration);
   atoms[i]->setVelocity(velocity);
   atoms[i]->setHighTimeDerivs(highTimeDerivs);
  }
}

// Method: gearCorrect
// Usage: gearCorrect(num, type, length, time, mass)
// ------------------------------------------------- 
// gearCorrect corrects the values of position,
// acceleration, velocity, third time derivative and
// forth time derivative predicted by the fifth order
// Gear method. It is called following calls to
// gearPredict and getForce.
 
double GearPC::gearCorrect(int num, double length, double delT)
{

  int i, j;
  double aNew;      // new accelerations
  double adjust;    // corrections to apply
  double mass, sumVSq;
  double c1, c2, c3, c4, c5;
  double cr, cv, cb, cc, cd;
  double *force, *acceleration, *position;
  double *velocity, **highTimeDerivs;

  c1 = delT;
  c2 = c1*delT/2;
  c3 = c2*delT/3;
  c4 = c3*delT/4;
  c5 = c4*delT/5;

  // Terms involving Gear's constants
  cr = G0*c2;
  cv = G1*c2/c1;
  cb = G3*c2/c3;
  cc = G4*c2/c4;
  cd = G5*c2/c5;
   
  sumVSq = 0.0;

  for(i = 0; i < num; i++)
  {
    mass = atoms[i]->getMass();
    acceleration = atoms[i]->getAcceleration();
    velocity = atoms[i]->getVelocity();
    highTimeDerivs = atoms[i]->getHighTimeDerivs();
 
   //determine corrections
   for(j = 0; j < 3; j++)
   {
     aNew = atoms[0]->f[3*i + j]/mass;
     adjust = aNew - acceleration[j];    //acceleration
     atoms[0]->r[3*i + j] += cr*adjust;
     atoms[0]->r[3*i + j] -= length * nearestInt(atoms[0]->r[3*i + j], length);
     velocity[j] += cv*adjust;
     acceleration[j] = aNew;
     highTimeDerivs[0][j] += cb*adjust;
     highTimeDerivs[1][j] += cc*adjust;
     highTimeDerivs[2][j] += cd*adjust;
   }
      
   atoms[i]->setAcceleration(acceleration);
   atoms[i]->setHighTimeDerivs(highTimeDerivs);

   // It is computationally convient to calculate the
   // kinetic energy too:
   sumVSq += mass*(velocity[0]*velocity[0] + velocity[1]*velocity[1]
                 + velocity[2]*velocity[2]);
   atoms[i]->setVelocity(velocity);     
  }

  return sumVSq/2;
}


// Constructor
// ------------


GearPC::GearPC(Atom **theAtoms)
:Integrator(theAtoms)
{
  //construct superclass
}


