// File:  atomMC.h
// -------------
// This file contains the definition of abstract class
// Atom, which defines all methods and data common to
// different atom types (eg derived classes), as well as
// foreshadowing methods of those derived classes through
// the use of abstract methods

#ifndef _atomMC_h_
#define _atomMC_h_

class Atom {
  protected:
    int type;         // type of atom
    double mass;      // mass of the atom
    double** rCutOff; // cut off distance for atom pairs
    double* position; // position
    double* trialPosition;

  public:
    int* kind;
    double* r; // public position vector array to aid GPU implementation
    double* rCutii;
    double* epsii;
    double* sigii;
    // Constructors
    // ------------
    // Used to construct the base class component of any
    // derived classes of Atom, and for the declaration of
    // an array of type Atom (for polymorphic method calls)

    Atom(int, double, int, int); // constructor
    Atom();                      // 2nd constructor for array declaration

    // Access (Set) Methods
    // --------------------
    void setType(int);         // assign type to atom
    void setrCutOff(double**); // assign cut off point for atom
    void setMass(double);      // assign mass to atom
    void setPosition(double*); // assign position of atom
    void setTrialPosition(double*);
    void reSetPosition();
    double** getrCutOff(); // return cut off point for atom
    double* getPosition(); // get position of atom
    double* getTrialPosition();
    double getMass(); // get mass of atom
    int getType();    // get type of atom

    // Abstract Virtual Methods
    // ------------------------
    // LJatom derived class methods
    // ----------------------------
    virtual void setSigma(double**) = 0;   // set sigma for atom pairs
    virtual void setEpsilon(double**) = 0; // set epsilon for atom pairs
    virtual double** getEpsilon() = 0;     // get epsilon for atom pairs
    virtual double** getSigma() = 0;       // get sigma for atom pairs
};

#endif
