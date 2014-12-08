/* state.cc    

Member functions of the class State 

Eric Jeckelmann, -- Uni Mainz -- October 2003

*/


/*
The class State describes a state of the hard-core particle system 
at a given time. 
*/

#include <cmath>
#include "state.h"
#include "setup.h"
#include <stdlib.h>


// initialization of the static class member 
Setup State::setup = Setup();


// return a pseudo-random number between -1/2 and 1/2
double randomize()
	{ return double(rand())/double(RAND_MAX)-0.5; }


/* 
Actual constructor  (construction of the initial state).
The box is divided in a regular array of cells. Each particle is placed alone in
a cell. The position within the cell is chosen randomly with a uniform probability. 
Initial velocities are randomly and uniformly distributed to make a state with 
vanishing total momentum and an average kinetic energy corresponding to the initial 
temperature.
*/
State::State(double t, const Setup & s)    
	{
	time = 0.0;
	temperature=t;
	setup = s;
	particles = vector<Particle>(setup.Np());

	//velocity bound for uniform velocity distribution
	double thermal_velocity=sqrt(3.0*temperature);

	//calculates cell size for initial position distribution
	int square_root_of_particle_number=int(sqrt(1.0*setup.Np()));
	if(square_root_of_particle_number*square_root_of_particle_number < setup.Np())
		square_root_of_particle_number++;
	double xCell_size=setup.xSize()/square_root_of_particle_number;
	double yCell_size=setup.ySize()/square_root_of_particle_number;

	double kinetic_energy=0.0;
	int i;
	int column = 0;
	int row = 0;
	for(i = 1; i < setup.Np(); i+=2)
		{
		double x = (column+0.5)*xCell_size+randomize()*(xCell_size-2*particle_radius);
		if(x < 1.0) x=1.0;
		if(x > setup.xSize()-1.0) x=setup.xSize()-1.0;
		double y = (row+0.5)*yCell_size+randomize()*(yCell_size-2*particle_radius);
		if(y < 1.0) y=1.0;
		if(y > setup.ySize()-1.0) y=setup.ySize()-1.0;
		column++;
		if(column == square_root_of_particle_number) 
			{
			column=0;
			row++;
			}
		double vx = randomize()*2*thermal_velocity;   
		double vy = randomize()*2*thermal_velocity;    
		particles[i-1] = Particle(x,y,vx,vy);
		x = (column+0.5)*xCell_size+randomize()*(xCell_size-2*particle_radius);
		if(x < 1.0) x=1.0;
		if(x > setup.xSize()-1.0) x=setup.xSize()-1.0;
		y = (row+0.5)*yCell_size+randomize()*(yCell_size-2*particle_radius);
		if(y < 1.0) y=1.0;
		if(y > setup.ySize()-1.0) y=setup.ySize()-1.0;
		column++;
		if(column == square_root_of_particle_number) 
			{
			column=0;
			row++;
			}
		particles[i] = Particle(x,y,-vx,-vy);      
		kinetic_energy+=vx*vx+vy*vy;
		}
	if(setup.Np()%2 == 1)
		{
		double x = (column+0.5)*xCell_size+randomize()*(xCell_size-2*particle_radius);
		if(x < 1.0) x=1.0;
		if(x > setup.xSize()-1.0) x=setup.xSize()-1.0;
		double y = (row+0.5)*yCell_size+randomize()*(yCell_size-2*particle_radius);
		if(y < 1.0) y=1.0;
		if(y > setup.ySize()-1.0) y=setup.ySize()-1.0;
		if(setup.Np() > 1)
			particles[setup.Np()-1] = Particle(x,y,0.0,0.0);      
		else
			{
			particles[setup.Np()-1] = Particle(x,y,sqrt(temperature),sqrt(temperature));
			kinetic_energy = temperature;
			}
		}
	double factor = temperature*setup.Np()/kinetic_energy;	
	for(i = 1; i <= setup.Np(); i++)
		{ particles[i-1].RescaleVelocity(factor); }  
	}

// Return the kinetic energy distribution
vector<int> State::statistic(int nbins, double resolution)
	{
	vector<int> result(nbins,0);
	for(int i=1; i<= setup.Np(); i++)
		{
		double energy = particles[i-1].kinetic_energy();
		int index = int(energy/resolution);
		if(index < nbins) result[index]++;
		}
	return result;
	}


// print out the current configuration
ostream & operator << (ostream & os, const State & state)
	{
	os << "Time = " << state.time << endl;
	for(int i=1; i<= state.Np(); i++)
		os << "Particle " <<  i << " : " << state.particles[i-1].Px() << " " << state.particles[i-1].Py()
		     << " " << state.particles[i-1].Vx() << " " << state.particles[i-1].Vy() << endl;
	return os;
	}

