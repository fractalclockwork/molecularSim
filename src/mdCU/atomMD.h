// File:  atomMD.h
// -------------
// This file contains the definition of abstract class
// Atom, which defines all methods and data common to
// different atom types (eg derived classes), as well as
// foreshadowing methods of those derived classes through
// the use of abstract methods

#ifndef _atomMD_h_
#define _atomMD_h_


class Atom
{
  private:
    int type;     		//type of atom
    double mass;     		//mass of the atom
    double *position;       	//position
    double *velocity;		//velocity
    double *acceleration;   	//acceleration
    double **highTimeDerivs;    //higher time derivatives
    double *force;   		//force
    double *forceOld;           //old force
  protected:
    double **rCutOff;           // cut off distance for atom pairs
  public:
    // Public variables to facilitate GPU implentation
    int *kind;
    double *r;
    double *f;
    double *rCutii;
    double *epsii;
    double *sigii;

    // Used to construct the base class component of any
    // derived classes of Atom, and for the declaration of
    // an array of type Atom (for polymorphic method calls)

    Atom(int, double, int, int, int);   //constructor
    Atom();                        //2nd constructor for array declaration

    //  Access (get and set) methods
    //  ----------------------------
 
    void setType(int);            	//assign type to atom
    void setrCutOff(double **);   	//assign cut off point for atom
    void setMass(double);         	//assign mass to atom
    void setAcceleration(double *); 	//assign acceleration to the atom
    void setPosition(double *);         //assign position of atom
    void setVelocity(double *);   	//assign velocity to atom
    void setHighTimeDerivs(double **);  //set derivatives of time
    void setForce(double *);      	//set the force of the atom
    void setForce(double,double,double);
    void setForce(int,double *);
    double **getrCutOff(); 		//return cut off point for atom
    double *getAcceleration();    	//get acceleration of atom
    double *getPosition();        	//get position of atom
    double *getVelocity();        	//get velocity of atom
    double **getHighTimeDerivs(); 	//get derivatives of time
    double *getForce();   		//get the force of the atom
    double *getForce(int);
    double getMass();             	//get mass of atom
    int getType();             		//get type of atom

    // Abstract Methods
    // ----------------
    // LJatom derived class methods
    // ----------------------------
    virtual void setSigma(double **) = 0;   //set sigma for atom pairs
    virtual void setEpsilon(double **) = 0; //set epsilon for atom pairs
    virtual double **getEpsilon() = 0;      //get epsilon for atom pairs
    virtual double **getSigma() = 0;        //get sigma for atom pairs
};

#endif


