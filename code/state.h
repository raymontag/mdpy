/* state.h   

Class State 

Eric Jeckelmann, -- Uni Mainz -- October 2003

*/

#ifndef _state_h_
#define _state_h_ 1

/*
Class describing the state of the hard-core particle system
at a given time in the molecular dynamics program. 
*/

#include <vector>
#include <cfloat>
using namespace std;
#include "particle.h"
#include "setup.h"


const double MAX_TIME = 0.01*DBL_MAX;
const double MIN_TIME = 10.0*DBL_EPSILON;


class State  
{
private:

	double time;     			// current time 
	double temperature;
	vector<Particle> particles; // array containing the particle positions and velocities 
	static Setup setup;   		// system setup

//make a copy 
	void copy(const State & other)
		{
		time = other.time;
		temperature = other.temperature;
		particles = other.particles; 
		}


public:

// Default constructor: null state   
	State()
		{ time = 0.0; }

// Actual constructor  (construction of the initial state)
	State(double t,  			// initial temperature
		  const Setup & s) ;

// Copy constructor
	State(const State & other)
		{ copy(other); }

// Copy operator / assignment operator
	State & operator=( State & other) 
		{
		if(this != &other) copy(other);
		return *this;
		}

// Return the system parameters

    double xSize() const
        { return setup.xSize(); }

    double ySize() const
        { return setup.ySize(); }

    int Np() const
        { return setup.Np(); }

// Return the kinetic energy distribution
	vector<int> statistic(int nbins, double resolution);

// print out the current configuration (positions and velocities)
	friend ostream & operator<<(ostream & os, const State & state);

// class used to compute the system dynamics
	friend class Dynamics;
};


#endif

