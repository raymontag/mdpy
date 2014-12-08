/* ising.cc    

Monte Carlo simulation for the two-dimensional Ising model
(demonstration program)

Eric Jeckelmann, Leibniz Universitaet Hannover 
Version 2, May 2009

-----

Definition of Ising model:

E({s_i}) = -J \sum_{<i,j>} s_i s_j - B \sum{i} s_i

s_i = +-1, ferromagnetic case J >=0 only.

----

Metropolis Monte Carlo simulation with deterministic sweeps
through the lattice and single spin flips only.

Units: J/kT (energy), B/kT (magnetic field), k=1 (Boltzmann constant)

NOTE: the random number generator seed is not initialized

*/

// uncomment to print out autocorrelation measurements
//#define _AUTOCORRELATIONS_

// uncomment to print sweep information
//#define _SWEEP_INFORMATION

// uncomment to print out spin dynamics   
//#define _SPIN_DYNAMICS_


#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
using namespace std;

int xsize;             // system length
int ysize;             // system width
int nspins;            // number of spins = system width * system length
double exchange;       // exchange coupling J >= 0
double field;          // magnetic field B
double temperature;    

vector <int> spins;    // spin configurations  
vector <int> previous_conf;

double energy;
double magnetization;


// return a pseudo-random number between 0 and 1
double randomize()
    { return double(rand())/double(RAND_MAX); }

//return the index of the spin with coordinates (ix,iy)
inline int index(int ix, int iy)
	{ return ix*ysize+iy; }

enum direction_t {Left,Up,Right,Down};

//return the index of the nearest-neighbors (periodic boundary conditions)
int neighbor(int ix, int iy, direction_t dir)
	{ 
	if(dir==Left)
		{
		ix--;
		if(ix == -1) ix=xsize-1;
		}
	else if(dir==Up)
		{
		iy++;
		if(iy == ysize) iy=0;
		}
	else if(dir==Right)
		{
		ix++;
		if(ix == xsize) ix=0;
		}
	else  // dir == Down
		{
		iy--;
		if(iy == -1) iy=ysize-1;
		}
	return index(ix,iy);
	}

// choose the initial spin configuration   
void initialization()
	{
	int order = -1;
	while(order != 0 && order != 1)
		{
		cerr << "Initial spin spins (0=random, 1=ordered)? ";
		cin >> order;
		cerr << endl;
		}
	spins.resize(nspins);
	if(order)
		{
		int spin = +1; 
		if(field < 0.0) spin =-1;
		cout << "Uniform initial spins with s=" << spin <<endl; 
		for(int i=0; i < nspins; i++) spins[i] = spin;   
		}
	else
		{
		cout << "Random initial spins" << endl; 
		for(int i=0; i < nspins; i++)
			{
			spins[i] = +1;
			if (randomize() < 0.5) spins[i] = -1;
			}
		}
	energy=0.0;
	for(int ix=0; ix < xsize; ix++)
	for(int iy=0; iy < ysize; iy++)
		{
		int i = index(ix,iy);
#ifdef _SPIN_DYNAMICS_
//		cout << "Spin: " << ix << " " << iy << " " << spins[i] << endl;
#endif
		energy-=exchange*spins[i]*spins[neighbor(ix,iy,Right)];
		energy-=exchange*spins[i]*spins[neighbor(ix,iy,Up)];
		energy-=field*spins[i];
		magnetization+=spins[i];
		}
	}


// main program
int main()
{

cerr << endl;
cerr << "Monte Carlo simulation for the two-dimensional Ising model" << endl;
cout << "Monte Carlo simulation for the two-dimensional Ising model" << endl;
cerr << endl;

// read in the system parameters  
xsize=ysize=0;
while(xsize < 3 && ysize < 3)
	{
	cerr << "System length and width? ";
	cin >> xsize >> ysize;
	cerr << endl;
	}
nspins=xsize*ysize;
cout << "System length and width = " << xsize << ", " << ysize << endl;

exchange=-1.0;
while(exchange < 0.0)
	{
	cerr << "Exchange coupling J, magnetic field B ? ";
	cin >> exchange >> field;
	cerr << endl;
	}

while(temperature <= 0.0)  
	{
	cerr << "Temperature (k=1)? ";
	cin >> temperature;
	cerr << endl;
	}
cout  << "Exchange coupling, magnetic field = " << exchange << ", " << field << endl;
cout << "Temperature = " << temperature << endl;
cout << endl;
exchange/=temperature;
field/=temperature;

// initialization
initialization();

cout << "Initial energy= " << energy*temperature << endl;
cout << "Initial magnetization= " << magnetization << endl;
cout << endl;


// possible probability ratio
vector<double> w_values(5); 
for(int i=0; i < 5; i++)
	w_values[i]=exp(2*exchange*(2*i-4)+2*field);

// read in the Monte Carlo simulation parameters
int sample_numbers=0;
int intermediate_sweeps=0;
int warmup_sweeps=-1;
while(warmup_sweeps < 0)
	{
	cerr << "Number of warmup sweeps? ";
	cin >> warmup_sweeps;
	cerr << endl;
	}
while(sample_numbers < 1)
	{
	cerr << "Number of samples? ";
	cin >> sample_numbers;
	cerr << endl;
	}
while(intermediate_sweeps < 1)
	{
	cerr << "Number of sweeps between samples? ";
	cin >> intermediate_sweeps;
	cerr << endl;
	}
int total_sweeps=max(warmup_sweeps,1)+intermediate_sweeps*(sample_numbers-1);
cout << warmup_sweeps << " warmup sweeps" << endl;
cout << sample_numbers << " samples" << endl;
cout << intermediate_sweeps << " sweeps between samples" << endl;
cout << "Total: " << total_sweeps << " sweeps" << endl;
cout << endl;


// statistical averages over many spin configurations computed with Monte Carlo
double average_flip_ratio = 0.0;
double average_autocorrel = 0.0;
double average_energy = 0.0;
double average_energy_square = 0.0;
double average_magnetization = 0.0;
double average_magnetization_square = 0.0;
int nsamples = 0;     // number of spin configurations used to compute the average

int ninter=intermediate_sweeps-1;
previous_conf = spins;
for(int ns=1; ns <= total_sweeps; ns++)
	{
	int nflips=0;    // number of spin flips in this sweep

// sweep through the lattice
	for(int ix=0; ix < xsize; ix++)
	for(int iy=0; iy < ysize; iy++)
		{
	// index of the spin that one tries to flip
		int i = index(ix,iy);
#ifdef _SPIN_DYNAMICS_
		if(nsamples > 0)
			cout << "Current: " << ix << " " << iy << " " << spins[i] << endl;
#endif
	// probability ratio between new and old spin configuration
	//	double ratio = exp(-2*field*spins[i]);   
		int nnspin = spins[neighbor(ix,iy,Left)];	//sum of n.n. spins	
		nnspin += spins[neighbor(ix,iy,Right)];	
		nnspin += spins[neighbor(ix,iy,Up)];	
		nnspin += spins[neighbor(ix,iy,Down)];	
		double ratio = w_values[(nnspin+4)/2];
		if(spins[i] == +1) ratio = 1.0/ratio;
		
	//Metropolis algorithm
		if(ratio < 1.0)
			{
			double z = randomize();
			if(ratio < z) continue;
			}

	// update for accepted spin flips
		nflips++;
		spins[i] *= -1;
		energy-=2*(exchange*nnspin+field)*spins[i];
		magnetization+=2*spins[i];
#ifdef _SPIN_DYNAMICS_
		if(nsamples > 0)
			cout << "Flip: " << ix << " " << iy << " " << spins[i] << endl;
#endif
		}

// printout information about this sweep
#ifdef _SWEEP_INFORMATION_
	cout << nflips << " spin flips (" << int(100*nflips/nspins) << "%) in sweep ";
	cout << ns << endl;
	cout << "Energy = " << energy*temperature << endl;
	cout << "Magnetization = " << magnetization << endl;
	cout << endl;
#endif

	if(ns < warmup_sweeps) continue;

	ninter++;
	if(ninter < intermediate_sweeps) continue;

	nsamples++;
	ninter=0;
	//cout << nsamples << " " << magnetization << endl;
	//cout << nsamples << " " << energy << endl;

// print out first spin configuration
#ifdef _SPIN_DYNAMICS_
	if(nsamples == 1)
		{
		for(int ix=0; ix < xsize; ix++)
		for(int iy=0; iy < ysize; iy++)
			{
			int i = index(ix,iy);
			cout << "Spin: " << ix << " " << iy << " " << spins[i] << endl;
			}
		}
#endif


// autocorrelation between spin configurations in successive samples
	double sum1 = 0;
	double sum2 = 0;
	double sum3 = 0;
	double sum4 = 0;
	double sum5 = 0;
	double autocorrel = 0.0;
	for(int j=0; j < nspins; j++)
		{
		sum1 += spins[j]*previous_conf[j];
		sum2 += spins[j]*spins[j];
		sum3 += previous_conf[j]*previous_conf[j];
		sum4 += spins[j];
		sum5 += previous_conf[j];
		}
	sum1 /= nspins;
	sum2 /= nspins;
	sum3 /= nspins;
	sum4 /= nspins;
	sum5 /= nspins;
	if(sum1 != sum4*sum5)
		autocorrel = (sum1-sum4*sum5)/sqrt(sum2-sum4*sum4)/sqrt(sum3-sum5*sum5);
	previous_conf = spins;

#ifdef _AUTOCORRELATIONS_
	cout << "Sample " << nsamples << endl;
	cout << "Autocorrelations = " << autocorrel << endl;
	if(sum2 == sum4*sum4 || sum3 == sum5*sum5)
		cout << "Warning: no fluctuation in spin configurations" << endl;
	cout << endl;
#else
	if(nsamples % 1000 == 0) cerr << "Sample " << nsamples << " done" << endl;
#endif

// calculate the MC statistical average
	average_flip_ratio += nflips; 
	average_autocorrel += abs(autocorrel);
	average_energy +=  energy;
	average_energy_square +=  energy*energy;
	average_magnetization +=  magnetization;
	average_magnetization_square +=  magnetization*magnetization;
	}


// print out the results
average_flip_ratio /= nsamples*nspins;
average_autocorrel /= nsamples;  
average_energy /= nsamples;
average_energy_square /= nsamples;
double energy_fluctuations = average_energy_square - average_energy * average_energy;
average_magnetization /=  nsamples;
average_magnetization_square /=  nsamples;
double magnetization_fluctuations = average_magnetization_square - average_magnetization * average_magnetization;
cout << "Number of samples = " << nsamples << endl;
cout << "Spin flip ratio = " << int(100*average_flip_ratio) << "%" << endl;
cout << "Autocorrelations = " << average_autocorrel << endl;
cout << "Energy per site = " << average_energy*temperature/nspins << endl;
cout << "Specific heat per site = " << energy_fluctuations/nspins << endl;
cout << "Magnetization per site = " << average_magnetization/nspins << endl;
cout << "Magnetic susceptibility per site = " << magnetization_fluctuations/temperature/nspins << endl;

if(exchange == 0.0)
	{
	double cosh = (exp(field)+exp(-field))/2.0;
	double sinh = (exp(field)-exp(-field))/2.0;
	double tanh = sinh/cosh;
	double sech = 1/cosh;
	cout << endl;
	cout << "Exact results: " << endl;
	cout << "Energy per site = " << -temperature*field*tanh << endl;
	cout << "Specific heat per site = " << field*field*sech*sech << endl; 
	cout << "Magnetization per site = " << tanh << endl;
	cout << "Magnetic susceptibility per site = " << sech*sech/temperature << endl;
	}

return 0;
}

