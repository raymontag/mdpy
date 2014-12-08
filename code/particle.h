/* particle.h

Class Particle

Eric Jeckelmann, -- Uni Mainz -- October 2003

*/

#ifndef _particle_h_
#define _particle_h_ 1

#include "setup.h"


// class describing a particle position and speed (velocity)

class Particle
{
private:

	double xPosition;
	double yPosition;
	double xVelocity;
	double yVelocity;

//make a copy 
	void copy(const Particle & other)
		{
		xPosition = other.xPosition;
		yPosition = other.yPosition;
		xVelocity = other.xVelocity;
		yVelocity = other.yVelocity;
		}

public:

// Default constructor: particle is immobile and out of system
	Particle()
		{ 
		xPosition = yPosition = -10*particle_radius; 
		xVelocity = yVelocity = 0.0;
		}

// Actual constructor
	Particle(double x, double y, double vx, double vy)
		{
		xPosition = x;
		yPosition = y;
		xVelocity = vx;
		yVelocity = vy;
		}

// Copy constructor
	Particle(const Particle & other)
		{ copy(other); }


// Copy operator / assignment operator
	Particle & operator=(const Particle & other) 
		{
		if(this != &other) copy(other);
		return *this;
		}

// Update position
	void Update(double time)
		{
		xPosition = xPosition + xVelocity * time;
		yPosition = yPosition + yVelocity * time;
		}

// Change particle speed
	void ChangeVelocity(double vx, double vy)
		{
		xVelocity = vx;
		yVelocity = vy;
		}

//Rescale particle speed
	void RescaleVelocity(double factor)
		{
		xVelocity *= factor;
		yVelocity *= factor;
		}

// Return the position and velocity vector components

	double Px() const
		{ return xPosition; }

	double Py() const
		{ return yPosition; }

	double Vx() const
		{ return xVelocity; }

	double Vy() const
		{ return yVelocity; }

// return the kinetic energy
	double kinetic_energy() const
		{ return 0.5*(xVelocity*xVelocity+yVelocity*yVelocity); }

};


#endif

