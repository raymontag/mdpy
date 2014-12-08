/* fastdynamics.cc

Member functions of the class Dynamics 

Eric Jeckelmann, -- Uni Mainz -- October 2003

*/

//#include <cfloat>
#include <cmath>
#include "fastdynamics.h"


/*
Motion of hard-code particles in two dimensions. The particles are confined in a 
rectangular box which encloses the points (x,y) with 0 < x < xLength and 0 < y < yLength.
The hard-core particle radius sets the unit length.
All collisions between hard-core particles are elastic (energy and momentum conservation).
Energy conservation and reflection law apply to collisions between a particle and one
of the box wall.
*/


//calculate the next wall hits and particle collisions
	void Dynamics::get_next_interactions()
		{
		int i;

		// calculate the next wall hit times
		wall_hit_times=vector<double>(Np());
		wall_indices=vector<int>(Np());
		first_wall_hit_time = MAX_TIME;
		for(i=1; i <= Np(); i++)
			{
			wall_hit_times[i-1]=time+next_wall_hit(particles[i-1],wall_indices[i-1]);
			if(first_wall_hit_time > wall_hit_times[i-1])
				{
				first_wall_hit_time = wall_hit_times[i-1];
				first_wall_hit_particle = i;
				}
			}

		// calculate the next particle collision times
		collision_times=vector<double>(Np());
		particle_indices=vector<int>(Np());
		first_collision_time = MAX_TIME;
		for(i=1; i <= Np(); i++)
			{
			next_collision(i); 
			if(first_collision_time > collision_times[i-1])
				{
				first_collision_time = collision_times[i-1];
				first_collision_particle = i;
				}
			}
		}

// return the time up to the next possible hit of the particle on a wall.
// wall_index designs the wall hit (1=bottom, 2=right, 3=top, 4=left)
	double Dynamics::next_wall_hit(const Particle & particle, int & wall_index)
		{
		double next_time_x=MAX_TIME;
		int wall_index_x=0;

	// particle hits the right wall
		if(particle.Vx() > 0.0) 
			{
			next_time_x = (xSize()-particle_radius-particle.Px())/particle.Vx();
			wall_index_x = 2;
			}
	// particle hits the left wall
		else if (particle.Vx() < 0.0) 
			{
			next_time_x = -(particle.Px()-particle_radius)/particle.Vx();
			wall_index_x = 4;
			}
		double next_time_y=MAX_TIME;
		int wall_index_y=0;
	// particle hits the top wall
		if(particle.Vy() > 0.0) 
			{
			next_time_y = (ySize()-particle_radius-particle.Py())/particle.Vy();
			wall_index_y = 3;
			}
	// particle hits the bottom wall
		else if (particle.Vy() < 0.0) 
			{
			next_time_y = -(particle.Py()-particle_radius)/particle.Vy();
			wall_index_y = 1;
			}

		double next_time = 0.0;
		if(next_time_x < next_time_y)
			{
			wall_index = wall_index_x;
			next_time = next_time_x;
			}
		else
			{
			wall_index = wall_index_y;
			next_time = next_time_y;
			}
		return next_time;
		}

// change the particle velocity after a collision with a wall.  
// Return the absolute momentum (i.e. velocity) change in the x-direction.
	double Dynamics::wall_hit(Particle & particle, int wall_index)
		{
		number_of_wall_hits++;
		double momentum=0.0;
		if(wall_index == 4)		// left wall
			{
			momentum = -2*particle.Vx();
			particle.ChangeVelocity(-particle.Vx(),particle.Vy());
			}
		if(wall_index == 2)		// right wall 
			{
			momentum = 2*particle.Vx();
			particle.ChangeVelocity(-particle.Vx(),particle.Vy());
			}
		if(wall_index == 1)		// lower wall 
			particle.ChangeVelocity(particle.Vx(),-particle.Vy());
		if(wall_index == 3)		// upper wall 
			particle.ChangeVelocity(particle.Vx(),-particle.Vy());
		return momentum;
		}

// calculate the next collision for a given particle
	void Dynamics::next_collision(int i)
		{
		double next_time = MAX_TIME;
		collision_times[i-1]= MAX_TIME;
		for(int j=1; j <= Np(); j++)
			{
			if(i==j) continue;
			double t = next_collision_time(particles[i-1],particles[j-1]);
			if(t < MIN_TIME) 
				cout << "Warning collision time = " << t 
					 << " for particles " << i << " and " << j << endl;
			if(t < next_time)
				{
				next_time=t;
				collision_times[i-1]=time+t;
				particle_indices[i-1]=j;
				}
			}
		}


// return the time up to the next possible collision between both particles  
	double Dynamics::next_collision_time(const Particle & particle1, const Particle & particle2)
		{
		double dvx=particle1.Vx()-particle2.Vx();
		double dvy=particle1.Vy()-particle2.Vy();
		double dpx=particle1.Px()-particle2.Px();
		double dpy=particle1.Py()-particle2.Py();
		double product=dvx*dpx+dvy*dpy;
		if(product >= 0.0) return MAX_TIME;
		double dv2=dvx*dvx+dvy*dvy;
		if(dv2 <= 0.0) return MAX_TIME;
		double dp2=dpx*dpx+dpy*dpy;
		double next_time=MAX_TIME;
		double determinant = product*product-dv2*(dp2-4.0*particle_radius*particle_radius);
		if(determinant > 0.0)
			next_time=(-product-sqrt(determinant))/dv2;
		else if(determinant == 0.0)
			{
			next_time = -product/dv2;
			if(next_time < MIN_TIME) next_time = MAX_TIME;  
			}
		return next_time;
		}

// change both particle speeds after a collision.
// Exchange the velocity components parallel to the vector between both particles.
	void Dynamics::collision(Particle & particle1, Particle & particle2)
		{
		number_of_collisions++;

		double dpx=particle1.Px()-particle2.Px();
		double dpy=particle1.Py()-particle2.Py();
		double factor1=(particle1.Vx()*dpx+particle1.Vy()*dpy)
		               /(4.0*particle_radius*particle_radius);
		double factor2=(particle2.Vx()*dpx+particle2.Vy()*dpy)
		               /(4.0*particle_radius*particle_radius);

		// components of parallel velocities in xy-basis
		double v1px=dpx*factor1;     
		double v1py=dpy*factor1;
		double v2px=dpx*factor2;     
		double v2py=dpy*factor2;

		// new velocity components in xy basis
		double v1x=v2px+particle1.Vx()-v1px;
		double v1y=v2py+particle1.Vy()-v1py;
		double v2x=v1px+particle2.Vx()-v2px;
		double v2y=v1py+particle2.Vy()-v2py;

		particle1.ChangeVelocity(v1x,v1y);
		particle2.ChangeVelocity(v2x,v2y);
		}


// let the state evolves until the next wall hit or particle collision
	void Dynamics::next_step(double movie_time_step)
		{
//update for all particles
		double time_step=first_collision_time-time; 
		if(first_wall_hit_time < first_collision_time)
			time_step=first_wall_hit_time-time;
		int i;
		if(time_step > movie_time_step) 
			{
			for(i=1; i <= Np(); i++)
				particles[i-1].Update(movie_time_step);
			time+=movie_time_step;
			return;
			}
		for(i=1; i <= Np(); i++)
			particles[i-1].Update(time_step);

//update for particle hitting the wall
		if(first_wall_hit_time < first_collision_time)
			{
			Particle & particle = particles[first_wall_hit_particle-1];
			momentum += wall_hit(particle,wall_indices[first_wall_hit_particle-1]);
			time = first_wall_hit_time;

			//update of wall hit information  
			wall_hit_times[first_wall_hit_particle-1]=
				time+next_wall_hit(particle,wall_indices[first_wall_hit_particle-1]);

			//update of collision information
			collision_times[first_wall_hit_particle-1] = MAX_TIME;
			for(int j=1; j <= Np(); j++)
				{
				if(first_wall_hit_particle == j) continue;
				double t = time+next_collision_time(particle, particles[j-1]);
				if(t < collision_times[first_wall_hit_particle-1])
					{
					collision_times[first_wall_hit_particle-1] = t;
					particle_indices[first_wall_hit_particle-1] = j;
					}
				if(particle_indices[j-1] == first_wall_hit_particle)
					next_collision(j);
				else if(t < collision_times[j-1])
					{
					collision_times[j-1] = t;
					particle_indices[j-1] = first_wall_hit_particle;
					}
				}
			}
//update for colliding particles 
		else
			{
			int index1 = first_collision_particle;    // index of first particle
			int index2 = particle_indices[index1-1];  // index of second particle
			Particle & particle1 = particles[index1-1];
			Particle & particle2 = particles[index2-1];

			collision(particle1,particle2);
			time = first_collision_time;

			//update of wall hit information  
			wall_hit_times[index1-1]=
				time+next_wall_hit(particle1,wall_indices[index1-1]);
			wall_hit_times[index2-1]=
				time+next_wall_hit(particle2,wall_indices[index2-1]);

			//update of collision information
			collision_times[index1-1] = MAX_TIME;
			collision_times[index2-1] = MAX_TIME;
			for(int j=1; j <= Np(); j++)
				{
				if(j == index1 || j == index2) continue;
				double t = time+next_collision_time(particle1, particles[j-1]);
				if(t < collision_times[index1-1])
					{
					collision_times[index1-1] = t;
					particle_indices[index1-1] = j;
					}
				if(particle_indices[j-1] == index1)
					next_collision(j);
				else if(t < collision_times[j-1])
					{
					collision_times[j-1] = t;
					particle_indices[j-1] = index1;
					}
				t = time+next_collision_time(particle2, particles[j-1]);
				if(t < collision_times[index2-1])
					{
					collision_times[index2-1] = t;
					particle_indices[index2-1] = j;
					}
				if(particle_indices[j-1] == index2)
					next_collision(j);
				else if(t < collision_times[j-1])
					{
					collision_times[j-1] = t;
					particle_indices[j-1] = index2;
					}
				}
			}

		first_wall_hit_time =MAX_TIME;
		first_collision_time =MAX_TIME;
		for(i=1; i <= Np(); i++)
			{
			if(first_wall_hit_time > wall_hit_times[i-1])
				{
				first_wall_hit_time = wall_hit_times[i-1];
				first_wall_hit_particle = i;
				}
			if(first_collision_time > collision_times[i-1])
				{
				first_collision_time = collision_times[i-1];
				first_collision_particle = i;
				}
			}
//		get_next_interactions();
		}

// let the state evolves between samplings
	void Dynamics::go_to_next_sampling(int make_movie, double movie_time_step)
		{
		init_counters();
		for(int i=1; number_of_wall_hits+number_of_collisions <= number_of_steps; i++)
			{
			next_step(movie_time_step);
			if(make_movie) cout << *this << endl;
			}
		get_next_interactions();   
		}


// print out information about the system (sampling)
	void Dynamics::printout()
		{
		cout << "Absolute time = " << time << endl;
		cout << "Time passed = " << time-initial_time << endl;
		cout << "Number of wall hits = " << number_of_wall_hits << endl;
		cout << "Number of particle collisions = " << number_of_collisions << endl;
		double interaction_time = (time-initial_time)/(number_of_collisions + number_of_wall_hits);
		cout << "Average time between interactions = " << interaction_time << endl;
		cerr << "Average time between interactions = " << interaction_time << endl;
		double pressure = momentum/(time-initial_time)/ySize()/2;
		cout << "Pressure = " << pressure << endl; 
//		cout << "Temperature = " << temperature << endl; 
		cout << "Ratio between temperature and pressure = " << temperature/pressure << endl;
		cout << endl;
		}




