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
  int i, k;
  double *position, *velocity, *acceleration, **highTimeDerivs;
  double c1, c2, c3, c4, c5;

  double sum = 0.0; 
  c1 = delT;
  c2 = c1*delT/2;
  c3 = c2*delT/3;
  c4 = c3*delT/4;
  c5 = c4*delT/5;
 
  for(i = 0; i < num; i++)
  {
    k = 3*i;

    velocity = atoms[i]->getVelocity();
    acceleration = atoms[i]->getAcceleration();
    highTimeDerivs = atoms[i]->getHighTimeDerivs();

    //new position vectors
    //position in the x-direction 
    atoms[0]->r[k]  += c1*velocity[0] + c2*acceleration[0]
		+  c3*highTimeDerivs[0][0]   //3rd derivative, x-direction
		+  c4*highTimeDerivs[1][0]   //4th derivative, x-direction 
		+  c5*highTimeDerivs[2][0];  //5th derivative, x-direction 
 
   //position in the y-direction 
    atoms[0]->r[k + 1]  += c1*velocity[1] + c2*acceleration[1]
	        +  c3*highTimeDerivs[0][1]   //3rd derivative, y-direction 
	        +  c4*highTimeDerivs[1][1]   //4th derivative, y-direction 
	        +  c5*highTimeDerivs[2][1];  //5th derivative, y-direction

    //position in the z-direction 
    atoms[0]->r[k + 2]  += c1*velocity[2] + c2*acceleration[2]
	        +  c3*highTimeDerivs[0][2]   //3rd derivative, z-direction
	        +  c4*highTimeDerivs[1][2]   //4th derivative, z-direction 
	        +  c5*highTimeDerivs[2][2];  //5th derivative, z-direction 

    atoms[0]->r[k]      -= length * nearestInt(atoms[0]->r[k], length);
    atoms[0]->r[k + 1]  -= length * nearestInt(atoms[0]->r[k + 1], length);
    atoms[0]->r[k + 2]  -= length * nearestInt(atoms[0]->r[k + 2], length);

    //new velocities 
    //x-componentt of velocity
    velocity[0] += c1*acceleration[0]
	        +  c2*highTimeDerivs[0][0] //3rd derivative, x-component
	        +  c3*highTimeDerivs[1][0] //4th derivative, x-component 
	        +  c4*highTimeDerivs[2][0];//5th derivative, x-component 

    //y-component of velocity
    velocity[1] += c1*acceleration[1]
	        +  c2*highTimeDerivs[0][1] //3rd Derivative, y-component 
 	        +  c3*highTimeDerivs[1][1] //4th Derivative, y-component
                +  c4*highTimeDerivs[2][1];//4th Derivative, y-component

    //z-component of  velocity
    velocity[2] +=  c1*acceleration[2]
                +   c2*highTimeDerivs[0][2]  //3rd derivative, z-component
                +   c3*highTimeDerivs[1][2]  //4th derivative, z-component 
                +   c4* highTimeDerivs[2][2]; //5th derivative, z-component

    //new accelerations
    //x-component of acceleration
    acceleration[0] += c1*highTimeDerivs[0][0] //3rd derivative, x-component 
		    +  c2*highTimeDerivs[1][0] //4th derivative, x-component
		    +  c3*highTimeDerivs[2][0];//5th derivative, x-component

    //y-component of acceleration
    acceleration[1] += c1*highTimeDerivs[0][1] //3rd derivative, y component
		    + c2*highTimeDerivs[1][1] //4th derivative, y component
                   + c3*highTimeDerivs[2][1];//5th derivative, y-component

    //z-component of acceleration
    acceleration[2] += c1*highTimeDerivs[0][2] //3rd derivative, z-component 
		    +  c2*highTimeDerivs[1][2] //4th derivative, z-component
		    +  c3* highTimeDerivs[2][2];//5th derivative, z-component 

    //new third derivatives
    //x-component of 3rd derivative
    highTimeDerivs[0][0] += c1*highTimeDerivs[1][0] //4th derivative, x-compt 
	                 +  c2* highTimeDerivs[2][0]; //5th derivative,x-compt


    //y-component of 3rd derivative
    highTimeDerivs[0][1] += c1*highTimeDerivs[1][1]  //4th der, y-compt
		         +  c2*highTimeDerivs[2][1]; //5th der, y-compt
 
    //z-component of 3rd derivative
    highTimeDerivs[0][2] += c1*highTimeDerivs[1][2]   //4th der, z-compt
		         +  c2*highTimeDerivs[2][2];  //5th der, z-compt
 

    //new fourth derivatives
    //x-component of 4th derivative
    highTimeDerivs[1][0] += c1*highTimeDerivs[2][0]; //5th der, x-compt
 
    //y-component of 4th derivative
    highTimeDerivs[1][1] += c1*highTimeDerivs[2][1]; //5th der, y-compt
 
    //z-component of 4th derivative 
    highTimeDerivs[1][2] += c1*highTimeDerivs[2][2]; //5th der, z-compt
     
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

  int i, k;
  double aNewX, aNewY, aNewZ;          // new accelerations
  double adjustX, adjustY, adjustZ;    // corrections to apply
  double mass, sumVSq;
  double c1, c2, c3, c4, c5;
  double cr, cv, cb, cc, cd;
  double *acceleration;
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
    k = 3*i;

    mass = atoms[i]->getMass();
    acceleration = atoms[i]->getAcceleration();
    velocity = atoms[i]->getVelocity();
    highTimeDerivs = atoms[i]->getHighTimeDerivs();
 
   //determine corrections
   aNewX = atoms[0]->f[k]/mass;      //x-force 
   aNewY = atoms[0]->f[k + 1]/mass;  //y-force
   aNewZ = atoms[0]->f[k + 2]/mass;  //z-force

   adjustX = aNewX - acceleration[0];    //x-acceleration
   adjustY = aNewY - acceleration[1];	 //y-acceleration
   adjustZ = aNewZ - acceleration[2];    //z-acceleration	
 
   //correct position vectors

   atoms[0]->r[k]     += cr*adjustX;
   atoms[0]->r[k + 1] += cr*adjustY;
   atoms[0]->r[k + 2] += cr*adjustZ;

   atoms[0]->r[k]     -= length * nearestInt(atoms[0]->r[k], length);
   atoms[0]->r[k + 1] -= length * nearestInt(atoms[0]->r[k + 1], length);
   atoms[0]->r[k + 2] -= length * nearestInt(atoms[0]->r[k + 2], length);
 
   //correct velocities
   velocity[0] += cv*adjustX;
   velocity[1] += cv*adjustY;
   velocity[2] += cv*adjustZ;
 
  //correct accelerations
  acceleration[0] = aNewX;
  acceleration[1] = aNewY;
  acceleration[2] = aNewZ;
 
  //correct third derivatives
  highTimeDerivs[0][0] += cb*adjustX;  
  highTimeDerivs[0][1] += cb*adjustY;  
  highTimeDerivs[0][2] += cb*adjustZ;  
 
  //correct fourth derivatives
  highTimeDerivs[1][0] += cc*adjustX; 
  highTimeDerivs[1][1] += cc*adjustY;  
  highTimeDerivs[1][2] += cc*adjustZ;  

  //correct fifth derivatives
  highTimeDerivs[2][0] += cd*adjustX;  
  highTimeDerivs[2][1] += cd*adjustY; 
  highTimeDerivs[2][2] += cd*adjustZ; 

  atoms[i]->setAcceleration(acceleration);
  atoms[i]->setVelocity(velocity);
  atoms[i]->setHighTimeDerivs(highTimeDerivs);

  // It is computationally convient to calculate the
  // kinetic energy too:
  sumVSq += mass*(velocity[0]*velocity[0] + velocity[1]*velocity[1] 
                 + velocity[2]*velocity[2]);
      
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


