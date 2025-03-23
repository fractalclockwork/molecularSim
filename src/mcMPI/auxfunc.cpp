// File: auxfuncMC.cpp
// -------------------_
// File containing auxilary methods common to many
// classes.

// This code was specifically developed to illustrate concepts in the accompanying book:
// R. J. Sadus, "Molecular Simulation of Fluids: Theory, Algorithms, Object-Orientation,
// and Parallel Computing," 2nd Ed. (Elsevier, Amsterdam, 2023). It can be used freely for
// any not for profit purpose or academic research application. The code has been validated,
// but it would be nonetheless prudent to test it further before publishing any results.
// Check the book's website for any subsequent updates.

#include <iostream>
using namespace std;

// Method: getMemory
// Usage: p = getMemory(num);
// -------------------------- 
// Assigns memory to 2-dimensional arrays.
 
double **getMemory(int num)
{
  int i;
  double **mat;
 
  mat = new double *[num];
 
  for(i = 0; i < num; i++)
    mat[i] = new double [num];
  
  return mat;
}
 
// Method: getMemory
// Usage: p = getMemory(xNum, yNum);
// --------------------------------- 
// Assigns memory to 2-dimensional arrays for higher time derivatives.
 
double **getMemory(int dim, int derivatives)
{
  int i;
  double **mat;

  mat = new double *[dim];
 
  for(i = 0; i < dim; i++)
    mat[i] = new double [derivatives];
  
  return mat;
}

// Method: nearestInt
// Usage: n = nearestInt(vector, boxL);
// ------------------------------------ 
// Return the next nearest integer resulting from
// the divison of vector by boxL.
 
double nearestInt(double vector, double boxL)
{
  int nInt;
  double div;
 
  div = vector/boxL;
 
  if(vector >= 0)
    nInt = (int)(div + 0.5);
  else
    nInt = -(int)(0.5 - div);
 
  return (double) nInt;
}

// Method: random
// Usage: n = random(idum);
// ------------------------ 
// Random initally takes any negative number (idum) and
// returns a random number on the interval 0 to 1. The
// algorithm is from Press et al.
 
double random(int *idum)
{

  //constants for random number generator
 
  const long int M1 = 259200;
  const int IA1 = 7141;
  const long int IC1 = 54773;
  const long int M2 = 134456;
  const int IA2 = 8121;
  const int IC2 = 28411;
  const long int M3 = 243000;
  const int IA3 = 4561;
  const long int IC3 = 51349;
 
  int j;
  static int iff = 0;
  static long int ix1, ix2, ix3;
  double RM1, RM2, temp;
  static double r[97];
 
  RM1 = 1.0/M1;
  RM2 = 1.0/M2;
 
  if(*idum < 0 || iff == 0)      //initialise on first call
  {  
    iff = 1;
    ix1 = (IC1 - (*idum)) % M1;   // seeding routines
    ix1 = (IA1 *ix1 + IC1) % M1;
    ix2 = ix1 % M2;
    ix1 = (IA1 * ix1 + IC1) % M1;
    ix3 = ix1 % M3;
 
    for (j = 0; j < 97; j++) 
    {
      ix1 = (IA1 * ix1 + IC1) % M1;
      ix2 = (IA2 * ix2 + IC2) % M2;
      r[j] = (ix1 + ix2 * RM2) * RM1;
     }
     *idum = 1;
  }
 
  //generate next number for each sequence
 
  ix1 = (IA1 * ix1 + IC1) % M1;
  ix2 = (IA2 * ix2 + IC2) % M2;
  ix3 = (IA3 * ix3 + IC3) % M3;
 
  //randomly sample r vector
 
  j = 0 + ((96 * ix3) / M3);
  if (j > 96 || j < 0)
     cout <<"Error in random number generator\n";
 
  temp = r[j];
 
  //replace r[j] with next value in the sequence
 
  r[j] = (ix1 + ix2 * RM2) * RM1;
 
  return temp;
}
