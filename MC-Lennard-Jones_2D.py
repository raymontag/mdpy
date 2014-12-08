# Monte Carlo simulation of a two-dimensional Lennard Jones gas
# Using periodic boundary conditions

# Eric Jeckelmann -- Leibniz Universitaet Hannover -- November 28, 2014


from visual import *
from visual.graph import *
import numpy as np
from random import uniform

debug=True

# visualization parameters
myrate=400
showEvolution=True     # show Monte Carlo dynamics
sweepsPerFrame=1       # update graphics after this number of sweeps

# Boltzmann constant in joules / kelvin 
#kB = 1.3806488e-23
# 1eV in joules
#eV = 1.60217657e-19

# Argon parameters
#mass=6.63e-26       # mass in kg
#sigma=3.822e-10     # position of LJ-potential minimum in meters
#radius=1.88e-10     # van der Waals radius in meters 
#sigma2=sigma*sigma
#TLJ=119.8           # LJ-potential minimum in kelvins
#epsilon=kB*TLJ      # LJ-potential minimum in joules 
mass=1.
sigma=1.
sigma2=sigma*sigma
radius=0.5
epsilon=1.
kB=1.

# system parameters
dimension=2          # box dimension (2=square, 3=cube)
N1=8
N=N1**dimension      # number of atoms
L=15.                # box edge size
V=L**dimension       # box volume
T=0.7                # temperature
idealP=N*kB*T/V
LJrangeCutOff=1.0/(2.5*sigma)**2   # inverse square of range of LJ potential


# Monte Carlo parameters
nTries=1000*N    # number of Metropolis steps
delta=0.6*sigma  # maximal step size (in units of sigma)
tauR=100*N       # relaxation time in number of MC steps (guess value)
tauC=100*N       # correlation time in number of MC steps (guess value)


# calculate the distances between all pairs in one direction
# for periodic boundary conditions
def Distances(r):
    d=r-r[:,newaxis]     # d[i,j] = r[i]-r[j]
    trim=greater(d,L/2)  # periodic boundary conditions  
    d=d-L*trim                  
    trim=less(d,-L/2)        
    d=d+L*trim
    return d

# Lennard-Jones potential energy with a finite interaction range
# and periodic boundary conditions
def PotentialEnergy(rx,ry,rz):
    dx=Distances(rx)
    dy=Distances(ry)
    dz=Distances(rz)
    d2=dx*dx+dy*dy+dz*dz + eye(N)
    invd2=sigma2*(1.00/d2 - eye(N))
    trim=greater(invd2,LJrangeCutOff)   # cut off potential at finite range
    invd2=invd2*trim
    invd2diff=invd2**6-2*invd2**3
    Vij=epsilon*invd2diff          # potential energy between each pair
    V=sum(Vij)/2
    return V

# calculate the energy change when one particle is moved
def EnergyChange(j,rxj,ryj,rzj,rx,ry,rz):
    # new position
    dx=rx-rxj
    trim=greater(dx,L/2)    
    dx=dx-L*trim                  
    trim=less(dx,-L/2)        
    dx=dx+L*trim
    dy=ry-ryj
    trim=greater(dy,L/2)    
    dy=dy-L*trim                  
    trim=less(dy,-L/2)        
    dy=dy+L*trim
    dz=rz-rzj
    trim=greater(dz,L/2)    
    dz=dz-L*trim                  
    trim=less(dz,-L/2)        
    dz=dz+L*trim
    d2=dx*dx+dy*dy+dz*dz
    d2[j]=1.
    invd2=sigma2/d2
    invd2[j]=0.
    trim=greater(invd2,LJrangeCutOff)  # cut off potential at finite range
    invd2=invd2*trim
    invd2diff=invd2**6-2*invd2**3
    Vij=epsilon*invd2diff
    newV=sum(Vij)
    # old position
    dx=rx-rx[j]
    trim=greater(dx,L/2)    
    dx=dx-L*trim                  
    trim=less(dx,-L/2)        
    dx=dx+L*trim
    dy=ry-ry[j]
    trim=greater(dy,L/2)    
    dy=dy-L*trim                  
    trim=less(dy,-L/2)        
    dy=dy+L*trim
    dz=rz-rz[j]
    trim=greater(dz,L/2)    
    dz=dz-L*trim                  
    trim=less(dz,-L/2)        
    dz=dz+L*trim
    d2=dx*dx+dy*dy+dz*dz
    d2[j]=1.
    invd2=sigma2/d2
    invd2[j]=0.
    trim=greater(invd2,LJrangeCutOff)   # cut off potential at finite range
    invd2=invd2*trim
    invd2diff=invd2**6-2*invd2**3
    Vij=epsilon*invd2diff
    oldV=sum(Vij)
    return newV-oldV

    
# calculate the pressure using the virial theorem
# for periodic boundary conditions and a finite interaction range
def Pressure(rx,ry,rz):
    dx=Distances(rx)
    dy=Distances(ry)
    dz=Distances(rz)
    d2=dx*dx+dy*dy+dz*dz + eye(N)
    invd2=sigma2*(1.00/d2 - eye(N))
    trim=greater(invd2,LJrangeCutOff)     # cut off forces at finite range
    invd2=invd2*trim
    invd2diff=invd2**7-invd2**4
    f=12.*epsilon/sigma2
    fx=f*dx*invd2diff  
    fy=f*dy*invd2diff
    fz=f*dz*invd2diff
    ax=sum(fx,axis=0)/mass 
    ay=sum(fy,axis=0)/mass
    az=sum(fz,axis=0)/mass
    p=idealP+(np.dot(ax,rx)+np.dot(ay,ry)+np.dot(az,rz))/(V*dimension)   
    return p

# initial configuration: cubic lattice
def init3D():
    d=L/N1
    rx=zeros(N)
    ry=zeros(N)
    rz=zeros(N)
    for iz in arange(N1): 
        for iy in arange(N1): 
            for ix in arange(N1):
                rx[ix+iy*N1+iz*N1*N1]=(0.5+ix)*d-L/2
                ry[ix+iy*N1+iz*N1*N1]=(0.5+iy)*d-L/2
                rz[ix+iy*N1+iz*N1*N1]=(0.5+iz)*d-L/2
    return rx,ry,rz

# initial configuration: square lattice
def init2D():
    d=L/N1
    rx=zeros(N)
    ry=zeros(N)
    rz=zeros(N) 
    for iy in arange(N1): 
        for ix in arange(N1):
            rx[ix+iy*N1]=(0.5+ix)*d-L/2
            ry[ix+iy*N1]=(0.5+iy)*d-L/2
    return rx,ry,rz


#### main program

# print out some information
print "Monte Carlo simulations for a Lennard-Jones potential"
print
print "Parameters"
print "Temperature= %4.2g" %T
print "Ideal gas pressure= %4.2g" %idealP
print "Number of tries= %6d" %nTries
print "Initial delta= %8.4f" %delta

# initialization
if dimension==2:
    rx,ry,rz=init2D()
else:
    print "ERROR: only dimension=2 is implemented"
    exit()
PE=PotentialEnergy(rx,ry,rz)
print "Initial potential energy per particle= %8.4f" %(PE/N)
    
# draw the initial configuration
if showEvolution:
    minXY=-L/2-radius
    maxXY=L/2+radius
    scene2 = display(title="Lennard-Jones MC",x=0,y=0,width=600,height=600,range=1.2*L/2,foreground=color.black,background=color.white,fov=0.01)
    scene2.select()
    # square box
    square = curve(pos=[(minXY,minXY),(minXY,maxXY),(maxXY,maxXY),(maxXY,minXY),(minXY,minXY)]) 
    Atoms=[]
    for i in arange(N):
        Atoms=Atoms+[sphere(pos=(rx[i],ry[i],rz[i]),radius=radius,color=color.blue)]
    Atoms[0].color=color.red

# plot for potential energy
graphics1=gdisplay(x=500,width=500,title="Potential energy",foreground=color.black,background=color.white) #energies
energyCurve1=gcurve(color=color.blue)
energyCurve2=gcurve(color=color.red)
# plot for pressure
graphics2=gdisplay(x=500,width=500,title="Pressure",foreground=color.black,background=color.white)
pressureCurve1=gcurve(color=color.blue)
pressureCurve2=gcurve(color=color.red)
# random walk plot
#if showEvolution:
#    track=curve(color=color.red)
#    track.append(pos=(rx[0],ry[0]))

# main loop over the Metropolis MC steps
nAccepted=0           # number of accepted Metropolis steps
averageEnergy=0.
averagePressure=0.
count=0               # number of measurements
for t in arange(nTries):
    rate(myrate)
    # generating a random step
    j=random.randint(N)                          # choose a particle randomly
    rxj=rx[j]+random.uniform(-delta/2,delta/2)   # random step in x-direction
    if rxj>L/2: rxj=rxj-L                        # periodic BCs
    if rxj<-L/2: rxj=rxj+L
    ryj=ry[j]+random.uniform(-delta/2,delta/2)   # random step in y-direction
    if ryj>L/2: ryj=ryj-L                        # periodic BCs
    if ryj<-L/2: ryj=ryj+L
    rzj=0.
    # Metropolis algorithm
    DE=EnergyChange(j,rxj,ryj,rzj,rx,ry,rz)
    w=exp(-DE/T)
    u=random.random()
    if u < w:           # step is accepted
        PE+=DE
        nAccepted+=1
        rx[j]=rxj
        ry[j]=ryj
        rz[j]=rzj
        if debug :
            en=PotentialEnergy(rx,ry,rz)
            if abs(en-PE) > 1e-12:
                print "ERROR in energy at step %6d, change= %8.4g, diff.= %8.4g" %(t,DE,en-PE)
    # measurement and visualization
    if mod(t,sweepsPerFrame*N)==0:
        # update the configuration graphics
        if showEvolution:
            for i in arange(N):   
                Atoms[i].pos=(rx[i],ry[i],rz[i])
#           track.append(pos=(rx[0],ry[0]))
        # calculate and plot energy and pressure
        if t > tauR:
            count=count+1
            en=PotentialEnergy(rx,ry,rz)
            energyCurve1.plot(pos=(t,en/N))
            averageEnergy+=en
            energyCurve2.plot(pos=(t,averageEnergy/(N*count)))
            if not debug and abs(en-PE) > 1e-12:
                print "ERROR in energy at step %6d, change= %8.4g, diff.= %8.4g" %(t,DE,en-PE)
            p=Pressure(rx,ry,rz)
            pressureCurve1.plot(pos=(t,p))
            averagePressure+=p
            pressureCurve2.plot(pos=(t,averagePressure/count))

# calculate averages and print out results
averageEnergy/=count
averagePressure/=count
print ""
print "Final results"
print "Number of accepted steps= %6d, acceptance ratio= %2.2g" % (nAccepted,(1.0*nAccepted)/nTries)
print "Potential energy per particle= %8.4f" % (averageEnergy/N)
print "Pressure= %8.4f" % (averagePressure)
