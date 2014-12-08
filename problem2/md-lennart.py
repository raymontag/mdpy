from visual import *
from random import uniform, randint
from math import sqrt, exp
from time import sleep
import numpy as np

# ANGSTROM = 10**(-10)
ANGSTROM = 1
ANGSTROM_SCALE = 10e10
#ANGSTROM_SCALE = 1

# UNIT = 1.6 * 10**-27
UNIT = 1
# UNIT_SCALE = (1./1.6)* 10**27
UNIT_SCALE = 1

global WALL_SIZE
RADIUS =  0.5 * ANGSTROM # [100 pm]
MAX_VELOCITY = .1
RATE = 10000

global T
N = 20
DT = 0.01
M = 1 * UNIT

kB = 1
# kB = 1
# EPSILON = 119.8 * kB # [K]
SIGMA =  1 * ANGSTROM # [100 pm]
EPSILON = 1
# SIGMA = 1

SYSTEM = "gas"
if SYSTEM == "gas":
    T = 1
    WALL_SIZE = 6
elif SYSTEM == "liquid":
    T = 0.0000001
    WALL_SIZE = 12
elif SYSTEM == "liquid2":
    T = 0.0000000001
    WALL_SIZE = 3.01
elif SYSTEM == "solid":
    T = 1e-60
    WALL_SIZE = 3.01


def dot(a, b):
    return sum(i*j for i,j in zip(a,b))

def vec_minus(a, b):
    return [i-j for i,j in zip(a,b)]

def vec_plus(a, b):
    return [i+j for i,j in zip(a,b)]

def norm(a):
    return sqrt(sum(i**2 for i in a))

def rescale(a, r):
    return [r*i for i in a]

def fit(x, a, b):
    return a*x+b

# Lennard-Jones-potential
def ljpot(r):
    #if r > RC:
    #    return 0
    return EPSILON * ((SIGMA / r)**12 - (SIGMA / r)**6)

class System(object):
    def __init__(self):
        self.balls = []

        # system = "gas"
        if SYSTEM == "gas":
            i = 0
            while i < N:
                x = uniform(0, WALL_SIZE)
                y = uniform(0, WALL_SIZE)
                z = uniform(0, WALL_SIZE)

                try:
                    for j in self.balls:
                        if norm(vec_minus([x, y, z], j.pos)) <= 2*RADIUS:
                            raise Exception
                except:
                    continue
                else:
                    self.balls.append(sphere(pos=(x,y,z), radius=RADIUS, color=color.red))
                    i += 1
        elif SYSTEM == "solid" or SYSTEM == "liquid2":
            x = 0.5
            while x <= 2.5:
                y = 0.5
                while y <= 2.5:
                    z = 0.5
                    while z <= 2.5:
                        self.balls.append(sphere(pos=(x,y,z), radius=RADIUS, color=color.red))
                        z += 1
                    y += 1
                x += 1
        else:
            x = 0.6
            while x <= 3.2:
                y = 0.6
                while y <= 3.2:
                    z = 0.6
                    while z <= 3.2:
                        self.balls.append(sphere(pos=(x,y,z), radius=RADIUS, color=color.red))
                        z += 1.2
                    y += 1.2
                x += 1.2


class Simulation(object):
    def __init__(self):
        self.scene = display(title='Lennard-Jones-MC',
                             x=0, y=0, width=300, height=300,
                             center=(WALL_SIZE/2,WALL_SIZE/2,WALL_SIZE/2),
                             background=(0,0,0), fullscreen = False)

        self.walls = []
        self.walls.append(box(pos=(WALL_SIZE/2,WALL_SIZE/2,0), length=WALL_SIZE,
                              height=WALL_SIZE, width=0.01, opacity=0.1))
        self.walls.append(box(pos=(WALL_SIZE/2,0,WALL_SIZE/2), length=WALL_SIZE,
                              height=0.01, width=WALL_SIZE, opacity=0.1))
        self.walls.append(box(pos=(0,WALL_SIZE/2,WALL_SIZE/2), length=0.01,
                              height=WALL_SIZE, width=WALL_SIZE, opacity=0.1))
        self.walls.append(box(pos=(WALL_SIZE/2,WALL_SIZE/2,WALL_SIZE), length=WALL_SIZE,
                              height=WALL_SIZE, width=0.01, opacity=0.1))
        self.walls.append(box(pos=(WALL_SIZE/2,WALL_SIZE,WALL_SIZE/2), length=WALL_SIZE,
                              height=0.01, width=WALL_SIZE, opacity=0.1))
        self.walls.append(box(pos=(WALL_SIZE,WALL_SIZE/2,WALL_SIZE/2), length=0.01,
                              height=WALL_SIZE, width=WALL_SIZE, opacity=0.1))

        self.system = System()
        self.old_system = []
        sleep(2)
        self.init_verlet()

    def run(self):
        while True:
            rate(RATE)
            self.do_verlet()

    def F(self, i):
        F = [0,0,0]
        j = 0
        while j < len(self.system.balls):
            if i == j:
                j += 1
                continue

            ri = np.copy(self.system.balls[i].pos)
            rj = np.copy(self.system.balls[j].pos)

            rijx = ri[0] - rj[0]
            if abs(rijx) > WALL_SIZE/2:
                rj[0] += WALL_SIZE * rijx / abs(rijx)
                rijx = ri[0] - rj[0]

            rijy = ri[1] - rj[1]
            if abs(rijy) > WALL_SIZE/2:
                rj[1] += WALL_SIZE * rijy / abs(rijy)
                rijy = ri[1] - rj[1]

            rijz = ri[2] - rj[2]
            if abs(rijz) > WALL_SIZE/2:
                rj[2] += WALL_SIZE * rijz / abs(rijz)
                rijz = ri[2] - rj[2]

            rij = (rijx, rijy, rijz)
                
            r = norm(rij)
            lj = ljpot(r)
            F = vec_plus(F, rescale(rij, 12 * lj / r**2))
            j += 1
        return F

    def init_verlet(self):
        init_v = []
        for i in self.system.balls:
            x = uniform(-MAX_VELOCITY, MAX_VELOCITY)
            y = uniform(-sqrt(MAX_VELOCITY**2 - x**2), sqrt(MAX_VELOCITY**2 - x**2))
            if randint(0, 1) == 0:
                z = -sqrt(MAX_VELOCITY**2 - x**2 - y**2)
            else: 
                z = sqrt(MAX_VELOCITY**2 - x**2 - y**2)
            init_v.append([x,y,z])

        i = 0
        new_system = []
        while i < len(self.system.balls):
            self.old_system.append(np.copy(self.system.balls[i].pos))
            new_system.append(vec_plus(vec_plus(self.system.balls[i].pos, rescale(init_v[i], DT)), rescale(self.F(i), DT**2 / (2*M))))
            i += 1

        i = 0
        while i < len(self.system.balls):
            self.system.balls[i].pos.x = new_system[i][0]
            self.system.balls[i].pos.y = new_system[i][1]
            self.system.balls[i].pos.z = new_system[i][2]
            i += 1

    def do_verlet(self):
        new_system = []
        i = 0
        while i < len(self.system.balls):
            new_system.append(vec_plus(vec_minus(rescale(self.system.balls[i].pos, 2), np.copy(self.old_system[i])), rescale(self.F(i), DT**2 / M)))
            i += 1

        # test = 0
        # old_system_copy = self.old_system[:]

        self.old_system = []
        i = 0
        while i < len(self.system.balls):
            if new_system[i][0] < 0:
                new_system[i][0] += WALL_SIZE
                self.system.balls[i].pos.x += WALL_SIZE
                # self.system.balls[i].pos.x %= WALL_SIZE
            elif new_system[i][0] > WALL_SIZE:
                new_system[i][0] -= WALL_SIZE
                self.system.balls[i].pos.x -= WALL_SIZE
                # self.system.balls[i].pos.x %= WALL_SIZE

            if new_system[i][1] < 0:
                new_system[i][1] += WALL_SIZE
                self.system.balls[i].pos.y += WALL_SIZE
                # self.system.balls[i].pos.y %= WALL_SIZE
            elif new_system[i][1] > WALL_SIZE:
                new_system[i][1] -= WALL_SIZE
                self.system.balls[i].pos.y -= WALL_SIZE
                # self.system.balls[i].pos.y %= WALL_SIZE

            if new_system[i][2] < 0:
                new_system[i][2] += WALL_SIZE
                self.system.balls[i].pos.z += WALL_SIZE
                # self.system.balls[i].pos.z %= WALL_SIZE
            elif new_system[i][2] > WALL_SIZE:
                new_system[i][2] -= WALL_SIZE
                self.system.balls[i].pos.z -= WALL_SIZE
                # self.system.balls[i].pos.z %= WALL_SIZE

            # Print velocity for test
            # test += (norm(vec_minus(new_system[i], old_system_copy[i])) / (2 * DT))**2

            self.old_system.append(np.copy(self.system.balls[i].pos))
            self.system.balls[i].pos.x = new_system[i][0]
            self.system.balls[i].pos.y = new_system[i][1]
            self.system.balls[i].pos.z = new_system[i][2]
            i += 1

        # print(test)
    
if __name__ == "__main__":
    simulation = Simulation()
    simulation.run()
