//  File: auxfuncMC.h
//  -----------------
//  Library interface for general purpose functions
//  that are not spectic to any class.

#ifndef _auxfuncMC_h
#define _auxfuncMC_h

int nearestInt(double vector, double boxL);
double **getMemory(int num);
double **getMemory(int xNum, int yNum);
double random(int *idum);
 
#endif
