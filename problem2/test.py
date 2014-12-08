from visual import *
from random import uniform, randint
from math import sqrt, exp
import numpy as np

WALL_SIZE = 50
RADIUS = 5
RATE = 10

N = 20
T = 50
DELTA = .1

# kB = 1.0e-23
EPSILON = 119.8 # * kB # [K * k_B]
SIGMA =  3.822 # [100 pm]

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
    return EPSILON * ((SIGMA / r)**12 - 2 * (SIGMA / r)**6)

class System(object):
    def __init__(self):
        self.balls = []

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


class Simulation(object):
    def __init__(self, delta):
        self.scene = display(title='Lennard-Jones-MC',
                             x=0, y=0, width=600, height=600,
                             center=(WALL_SIZE/2,WALL_SIZE/2,WALL_SIZE/2),
                             background=(0,0,0))

        self.walls = []
        self.walls.append(box(pos=(WALL_SIZE/2,WALL_SIZE/2,-1), length=WALL_SIZE,
                              height=WALL_SIZE, width=2, opacity=0.1))
        self.walls.append(box(pos=(WALL_SIZE/2,-1,WALL_SIZE/2), length=WALL_SIZE,
                              height=2, width=WALL_SIZE, opacity=0.1))
        self.walls.append(box(pos=(-1,WALL_SIZE/2,WALL_SIZE/2), length=2,
                              height=WALL_SIZE, width=WALL_SIZE, opacity=0.1))
        self.walls.append(box(pos=(WALL_SIZE/2,WALL_SIZE/2,WALL_SIZE+1), length=WALL_SIZE,
                              height=WALL_SIZE, width=2, opacity=0.1))
        self.walls.append(box(pos=(WALL_SIZE/2,WALL_SIZE+1,WALL_SIZE/2), length=WALL_SIZE,
                              height=2, width=WALL_SIZE, opacity=0.1))
        self.walls.append(box(pos=(WALL_SIZE+1,WALL_SIZE/2,WALL_SIZE/2), length=2,
                              height=WALL_SIZE, width=WALL_SIZE, opacity=0.1))


        self.system = System()
        self.w = self.init_w()
        self.delta = delta

    def init_w(self):
        w = 0
        for j in range(N):
            for i in range(j):
                xi = self.system.balls[i].pos[0]
                yi = self.system.balls[i].pos[1]
                zi = self.system.balls[i].pos[2]
                xj = self.system.balls[j].pos[0]
                yj = self.system.balls[j].pos[1]
                zj = self.system.balls[j].pos[2]

                if abs(xi - xj) < WALL_SIZE/2:
                    dx = abs(xi - xj)
                else:
                    dx = WALL_SIZE - abs(xi - xj)

                if abs(yi - yj) < WALL_SIZE/2:
                    dy = abs(yi - yj)
                else:
                    dy = WALL_SIZE - abs(yi - yj)

                if abs(zi - zj) < WALL_SIZE/2:
                    dz = abs(zi - zj)
                else:
                    dz = WALL_SIZE - abs(zi - zj)

                r = (dx, dy, dz)
                w += ljpot(norm(r))
        return w

    def calc_new_w(self, coords, i):
        dw = 0
        for j in range(N):
            if i == j:
                continue
            xi_new = coords[i][0]
            yi_new = coords[i][1]
            zi_new = coords[i][2]
            xi_old = self.system.balls[i].pos[0]
            yi_old = self.system.balls[i].pos[1]
            zi_old = self.system.balls[i].pos[2]
            xj = coords[j][0]
            yj = coords[j][1]
            zj = coords[j][2]

            if abs(xi_new - xj) < WALL_SIZE/2:
                dx = abs(xi_new - xj)
            else:
                dx = WALL_SIZE - abs(xi_new - xj)

            if abs(yi_new - yj) < WALL_SIZE/2:
                dy = abs(yi_new - yj)
            else:
                dy = WALL_SIZE - abs(yi_new - yj)

            if abs(zi_new - zj) < WALL_SIZE/2:
                dz = abs(zi_new - zj)
            else:
                dz = WALL_SIZE - abs(zi_new - zj)

            r = (dx, dy, dz)
            dw += ljpot(norm(r))

            if abs(xi_old - xj) < WALL_SIZE/2:
                dx = abs(xi_old - xj)
            else:
                dx = WALL_SIZE - abs(xi_old - xj)

            if abs(yi_old - yj) < WALL_SIZE/2:
                dy = abs(yi_old - yj)
            else:
                dy = WALL_SIZE - abs(yi_old - yj)

            if abs(zi_old - zj) < WALL_SIZE/2:
                dz = abs(zi_old - zj)
            else:
                dz = WALL_SIZE - abs(zi_old - zj)

            r = (dx, dy, dz)
            dw -= ljpot(norm(r))

        return self.w + dw

    def run(self):
        while True:
            rate(RATE)
            self.do_sweep()

    def do_metropolis(self):
        i = randint(0, N - 1)
        dx = uniform(-self.delta, self.delta)
        dy = uniform(-self.delta, self.delta)
        dz = uniform(-self.delta, self.delta)

        new_balls = [j.pos for j in self.system.balls]
        new_balls[i] = vec_plus(new_balls[i], (dx,dy,dz))
        new_w = self.calc_new_w(new_balls, i)
        print("NEW_W: "+str(new_w))
        print("OLD_W: "+str(self.w))
        print("DIFF: "+str(new_w-self.w))

        p = np.exp((self.w - new_w)/T)
        r = uniform(0,1)

        print("P: "+str(p))

        if p > r:
            self.w = new_w

            if new_balls[i][0] < 0:
                new_balls[i][0] = WALL_SIZE
            elif new_balls[i][0] > WALL_SIZE:
                new_balls[i][0] = 0

            if new_balls[i][1] < 0:
                new_balls[i][1] = WALL_SIZE
            elif new_balls[i][1] > WALL_SIZE:
                new_balls[i][1] = 0

            if new_balls[i][2] < 0:
                new_balls[i][2] = WALL_SIZE
            elif new_balls[i][2] > WALL_SIZE:
                new_balls[i][2] = 0

            self.system.balls[i].pos = new_balls[i]

            return True
        return False


    def do_sweep(self):
        successful_steps = 0.
        for i in range(N):
            if self.do_metropolis() is True:
                successful_steps += 1.

        success_rate = successful_steps/N
        if success_rate < 0.4:
            self.delta += 0.0001
        elif success_rate > 0.6:
            self.delta -= 0.0001
        print(success_rate)
                
    
if __name__ == "__main__":
    simulation = Simulation(DELTA)
    simulation.run()
