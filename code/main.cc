/* main.cc    

Molecular dynamics program for hard-core particles in two dimensions.

Eric Jeckelmann -- University of Mainz, October 2003
				-- Leibniz Universitaet Hannover, July 2011


The system consists of hard-core particles in two dimensions.
The particles are confined in a rectangular box.  
All collisions between particles and between a particle and the box walls
are elastic (=> conservation of energy and momentum).

The hard-core particle radius set the unit length.
Mass and Boltzmann constant are set to 1.
*/


#include <iostream>
using namespace std;

#include "dynamics.h"


int main()
{

cerr << endl;
cerr << "Hard-core particle molecular dynamics" << endl;
cout << "Hard-core particle molecular dynamics" << endl;
cerr << endl;

Setup setup;
setup.readin();

// read in the initial system temperature ( > 0)
double initial_temperature = 0.0;
while(initial_temperature <= 0.0)  
	{
	cerr << "Initial system temperature? ";
	cin >> initial_temperature;
	cerr << endl;
	}
cout << "Initial temperature = " << initial_temperature << endl;
cout << endl;

// read in the movie flag and time step ( > 0)
int make_movie = 0;
double movie_time_step = 0.0;
cerr << "Save configurations for movie (type 1 for yes)? ";
cin >> make_movie;
cerr << endl;
if(make_movie)
	{
	cout << "Saving configurations for movie" << endl;
	while(movie_time_step <= 0.0)  
		{
		cerr << "Movie time step? ";
		cin >> movie_time_step;
		cerr << endl;
		}
	cout << "Movie time step = " << movie_time_step << endl;
	cout << endl;
	}
else
	{
	cout << "Configurations not saved" << endl;
	movie_time_step = MAX_TIME;
	}


State state(initial_temperature,setup);
Dynamics system(state);
if(make_movie) 
	{
	cout << state << endl;
	system.go_to_next_sampling(make_movie,movie_time_step);
	system.printout();
	}
else
	{
	vector<int> distribution = state.statistic(100,initial_temperature/10);
	cout << "Kinetic energy distribution (resolution = " << initial_temperature/10 << ")" << endl;
	for(int i=0; i<100; i++)
		cout << "Bin " << i << " " << distribution[i] << endl;
	cout <<endl;
	for(int i=1; i<=10; i++)
		{
		system.go_to_next_sampling(make_movie,movie_time_step);
		cout << "Sample " << i << endl;
		system.printout();
		}
	distribution = system.state().statistic(100,initial_temperature/10);
	cout << "Kinetic energy distribution (resolution = " << initial_temperature/10 << ")" << endl;
	for(int i=0; i<100; i++)
		cout << "Bin " << i << " " << distribution[i] << endl;
	cout <<endl;
	}

return 0;
}

