/* fastdynamics.h 

Class Dynamics for describing the motion of hard-core particles.

Eric Jeckelmann, -- Uni Mainz -- October 2003

*/

#ifndef _dynamics_h_
#define _dynamics_h_ 1

#include "particle.h"
#include "state.h"


/*
Motion of hard-code particles in two dimensions. The particles are confined in a 
rectangular box which encloses the points (x,y) with 0 < x < xLength and 0 < y < yLength.
The hard-core particle radius sets the unit length.
All collisions between hard-core particles are elastic (energy and momentum conservation).
Energy conservation and reflection law apply to collisions between a particle and one
of the box wall.
*/


class Dynamics : private State
{
private:

	double initial_time;				// initial time of this simulation

// Information about next wall hits 
	vector<double> wall_hit_times;      // for each particle
	vector<int> wall_indices;  	    	// wall hit (1=bottom,2=right,3=top,4=left) 
	double first_wall_hit_time;         // time of first wall hit
	int first_wall_hit_particle;        // index of first particle hitting the wall

// Information about next particle collisions 
	vector<double> collision_times;     // for each particle
	vector<int> particle_indices;  	   	// index of collision partner
	double first_collision_time;        // time of first collision 
	int first_collision_particle;       // index of one particle involved in the first collision 

// Hit wall and collision counts
	long number_of_wall_hits;			
	long number_of_collisions;			

// Total momentum change at the wall (used to calculate the pressure
	double momentum;

// Number of time steps (hits and collisions) between measurements of the system
// properties
	int number_of_steps;


//private function members

// initialize counters, etc ...
	void init_counters()
		{
		momentum=0.0;
		number_of_wall_hits=0;
		number_of_collisions=0;
		initial_time = time;
		}

// determine the next wall hits and particle collisions
	void get_next_interactions();

// return the time up to the next possible hit of the particle on a wall.
// wall_index designs the wall hit (1=bottom, 2=right, 3=top, 4=left)
	double next_wall_hit(const Particle & particle, int & wall_index);

// change the particle velocity after a collision with a wall. 
// Return the (absolute) momentum (i.e. velocity) change in the x-direction.
	double wall_hit(Particle & particle, int wall_index);

// calculate the next collision for a given particle  
	void next_collision(int i);

// return the time up to the next possible collision between both particles  
	double next_collision_time(const Particle & particle1, const Particle & particle2);

// change both particle speeds after a collision.
// Exchange the velocity components parallel to the vector between both particles.
	void collision(Particle & particle1, Particle & particle2);

// let the state evolves until the next wall hit or particle collision
	void next_step(double movie_time_step);


//hidden automatic members

//default constructor
	Dynamics() { abort(); };

//copy constructor
	Dynamics(const Dynamics & other) { abort(); };

//copy/assignment operator
	Dynamics & operator=(const Dynamics & other) { abort(); return *this; };


public:

//actual constructor (initialization)
	Dynamics(const State & state): State(state)
		{
		number_of_steps = max(30,3*Np());
		init_counters();
		get_next_interactions();
		}

// let the system evolves between samplings
	void go_to_next_sampling(int make_movie, double movie_time_step);

// return a reference to the current system state
	State & state()
		{ return *this; }

// print out information about the system (sampling)
	void printout();

};

#endif

