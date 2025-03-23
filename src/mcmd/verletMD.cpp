// File: verletMD.cpp
// Class:Verlet 
 
// This code was specifically developed to illustrate concepts in the accompanying book:
// R. J. Sadus, "Molecular Simulation of Fluids: Theory, Algorithms, Object-Orientation,
// and Parallel Computing," 2nd Ed. (Elsevier, Amsterdam, 2023). It can be used freely for
// any not for profit purpose or academic research application. The code has been validated,
// but it would be nonetheless prudent to test it further before publishing any results. 
// Check the book's website for any subsequent updates.

#include "verletMD.h"
#include "auxfunc.h"
#include <iostream>

using namespace std;

// Method: velVerletP1 
// Usage: velVerletP1(num, length);
// -------------------------------- 
// The function is the first part of the velocity Verlet  algorithm
// to predict new coordinates, and  velocities.

void Verlet::velVerletP1(int num, double length, double delT)
{
  int  i, old = 1;
  double *position, *velocity, *force;
  double mass, dt2, dtsq2, ax, ay, az;

  dt2  = delT/2.0;
  dtsq2  = delT*dt2;
 
  for(i = 0; i < num; i++)
  {
    mass = atoms[i]->getMass(); 
    position = atoms[i]->getPosition();
    velocity = atoms[i]->getVelocity();
    force  = atoms[i]->getAcceleration();

    //determine accelerations
    ax = force[0]/mass;
    ay = force[1]/mass;
    az = force[2]/mass;
 
    //new position vectors
    position[0] += delT*velocity[0] + dtsq2*ax;
    position[1] += delT*velocity[1] + dtsq2*ay; 
    position[2] += delT*velocity[2] + dtsq2*az; 

    position[0] -= length * nearestInt(position[0], length);
    position[1] -= length * nearestInt(position[1], length);
    position[2] -= length * nearestInt(position[2], length);

    //Store old forces for the second part of velocity-Verlet
    atoms[i]->setForce(1,force);

    //Update positions;
    atoms[i]->setPosition(position);
  }
}

// Method: velVerletP2 
// Usage:  x = velVerlet2(num, type, length, time, mass)
// ------------------------------------------------- 
//  Second part of the velocity-Velert integrator that
//  ajdusts the velocities following updated forces
//  folling a second call to setForces. It returns the
//  kinetic energy calculated using the updated velocities..
 
double Verlet::velVerletP2(int num, double length, double delT)
{

  int i;
  double mass, dt2, sumVSq;
  double *force, *forceOld;
  double *velocity;


  dt2  = delT/2.0;
  sumVSq = 0.0;

  for(i = 0; i < num; i++)
  {
    mass = atoms[i]->getMass();
    force = atoms[i]->getForce();
    forceOld = atoms[i]->getForce(1);
    velocity = atoms[i]->getVelocity();
 
   //adjust velocities
 
    velocity[0] += (force[0] + forceOld[0])*dt2/mass;
    velocity[1] += (force[1] + forceOld[1])*dt2/mass;
    velocity[2] += (force[2] + forceOld[2])*dt2/mass;

//    cout << force[0] << " " << forceOld[0] << endl;
   
    atoms[i]->setVelocity(velocity);

    // It is computationally convient to calculate the
    // kinetic energy too:
  
    sumVSq += mass*(velocity[0]*velocity[0] + velocity[1]*velocity[1] 
                 + velocity[2]*velocity[2]);
      
  }

  return 0.5*sumVSq;;
}


// Constructor
// ------------

Verlet::Verlet(Atom **theAtoms)
:Integrator(theAtoms)
{
  //construct superclass
}


