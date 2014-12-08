/* setup.h   

Class Setup 

Eric Jeckelmann, -- Uni Mainz -- October 2003

*/

#ifndef _setup_h_
#define _setup_h_ 1

#include <iostream>
using namespace std;

/*
Class describing the system setup for the molecular dynamics program. 

The system consists of hard-core particles in  two dimensions.
The particles are confined in a rectangular box which encloses the points 
(x,y) with 0 < x < xLength and 0 < y < yLength.
*/

// the hard-core particle radius set the unit length
const double particle_radius = 1.0;

// minimal free space around a particle (fraction)
const double free_space = 0.05;

// minimal volume per particle = [2 * particle_radius * (1 + free volume)]^2 
const double particle_volume = 4.41;   


class Setup
{
private:

	double xLength;     // box size 
	double yLength;
	int Nparticles;     // number of particles  

//make a copy 
	void copy(const Setup & other)
		{
		xLength = other.xLength;
		yLength = other.yLength;
		Nparticles = other.Nparticles; 
		}

public:

// Default constructor: null system 
	Setup()
		{ 
		xLength = yLength = 0.0; 
		Nparticles = 0;   
		}

// Copy constructor
	Setup(const Setup & other)
		{ copy(other); }


// Copy operator / assignment operator
	Setup & operator=(const Setup & other) 
		{
		if(this != &other) copy(other);
		return *this;
		}

// Read in parameters (initialization)  
	void readin();

// Return the parameters   

	double xSize() const
		{ return xLength; }

	double ySize() const
		{ return yLength; }

	int Np() const
		{ return Nparticles; }

};


#endif

