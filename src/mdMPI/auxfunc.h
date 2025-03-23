// File: auxfuncMD.h
//  ---------------
// Library interface for general purpose functions
// that are not spectic to any class.

#ifndef _auxfuncMD_h
#define _auxfuncMD_h

double** getMemory(int num);
double** getMemory(int xNum, int yNum);
int nearestInt(double vector, double boxL);
double random(int* idum);
void freeMemory(double**);
void freeMem(double**, int);

#endif
