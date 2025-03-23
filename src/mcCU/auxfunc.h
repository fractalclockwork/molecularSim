//  File: auxfuncMC.h
//  -----------------
//  Library interface for general purpose functions
//  that are not spectic to any class.

#ifndef _auxfuncMC_h
#define _auxfuncMC_h

double **getMemory(int num);
double **getMemory(int xNum, int yNum);
double nearestInt(double vector, double boxL);
double random(int *idum);
 
#endif
