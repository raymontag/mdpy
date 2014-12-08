import math 

class test(object):
    def init_wall_collides(self):
        self.wall_collides = []
        for i in range(len(self.coordinates)):
            x = self.coordinates[i][0]
            vx = self.velocities[i][0]

            if vx > 0:
                tx = (SCREEN_SIZE - (RADIUS + x)) / vx
            elif vx < 0:
                tx = (RADIUS - x) / vx
            else:
                tx = INFINITY

            y = self.coordinates[i][1]
            vy = self.velocities[i][1]

            if vy > 0:
                ty = (SCREEN_SIZE - (RADIUS + y)) / vy
            elif vy < 0:
                ty = (RADIUS - y) / vy
            else:
                ty = INFINITY

            self.wall_collides.append([tx, ty])
        print(self.wall_collides)

    def init_particle_collides(self):
        self.coordinates = [[358.19296276283137, 65.03740457221139], [323.78637897931173, 146.61285232729853]]
        self.velocities = [[-0.6068027428719919, -0.06979430069741088], [-0.7887001600249275, 0.069794300697411]]
        self.particle_collides = []

        for i in range(len(self.coordinates)):
            times = []

            r1 = self.coordinates[i]
            v1 = self.velocities[i]
            for j in range(i):
                if i == j:
                    times.append(INFINITY)
                    continue
                times.append(self.particle_collides[i-1][j])
            for j in range(i, len(self.coordinates)):
                if i == j:
                    times.append(INFINITY)
                    continue
                r2 = self.coordinates[j]
                v2 = self.velocities[j]
                r12 = vec_minus(r1, r2)
                v12 = vec_minus(v1, v2)
        
                d = (dot(r12, v12)**2 - dot(v12, v12) * (dot(r12, r12) - 4 * RADIUS**2))
                if d > 0:
                    t = (-dot(r12, v12)-math.sqrt(d))/dot(v12, v12)
                    if t > 0:
                        times.append(t)
                    else:
                        times.append(INFINITY)
                else:
                    times.append(INFINITY)

            self.particle_collides.append(times)
        print(self.particle_collides)
        
N = 2
RADIUS = 40
MAX_VELOCITY = 1
DT = 3
X_DIRECTION = 0
Y_DIRECTION = 1
WALL_COLLISION = -1

SCREEN_SIZE = 400
CPU_TICKS = 100
INFINITY = 99999999999

def dot(a, b):
    return sum(i*j for i,j in zip(a,b))

def vec_minus(a, b):
    return [i-j for i,j in zip(a,b)]

def vec_plus(a, b):
    return [i+j for i,j in zip(a,b)]

def norm(a):
    return math.sqrt(sum(i**2 for i in a))

def rescale(a, r):
    return [r*i for i in a]
                                                                               
a = test()
a.init_particle_collides()
a.init_wall_collides()
