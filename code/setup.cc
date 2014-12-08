/* setup.cc

Member of functions of the class Setup 

Eric Jeckelmann, -- Uni Mainz -- October 2003

*/

#include <iostream>
using namespace std;
#include "setup.h"


// Read in system parameters (initialization): box size and number of particles. 
// The box must be large enough for at least one particle.
// The number of particles must be between 1 and the maximal number of particles
// that the box can contain (=  box volume / particle volume).

void Setup::readin()
	{
	double minimal_length = 2*particle_radius*(1+free_space);
	while(xLength < minimal_length || yLength < minimal_length)
		{
		cerr << "System size Lx,Ly (minimum " << minimal_length << ")? ";
		cin >> xLength >> yLength;
		}
	int max_particles = int(xLength * yLength / particle_volume);
	while(Nparticles < 1 || Nparticles > max_particles)
		{
		cerr << "Number of particles (1 to " << max_particles << ")? ";
		cin >> Nparticles;
		}
	cerr << endl;
	cout << "System size = " << xLength << " x " << yLength << endl;
	cout << "Number of particles = " << Nparticles << endl;
	double density = Nparticles / (xLength * yLength);
	cout << "Density = " << density << endl;
	cout << "Volume per particle = " << 1/density << endl;
	cout << "Reduced volume per particle = " 
	     << 1/density-3.14159*particle_radius*particle_radius << endl;
	}

