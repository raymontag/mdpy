import pygame
import random
import math
import pylab as P
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit

WHITE = (255, 255, 255)
BLUE =  (  0,   0, 255)
GREEN = (  0, 255,   0)
RED =   (255,   0,   0)

N = 100
RADIUS = 5
MAX_VELOCITY = 2
DT = 1
X_DIRECTION = 0
Y_DIRECTION = 1
WALL_COLLISION = -1
EQUILIBRIUM = True

SCREEN_SIZE = 800
CPU_TICKS = 500
INFINITY = 99999999999
WRITE_PDF = False
PRESSURE = False
PRESSURE_LENGTH = 1000

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

def fit(x, a, b):
    return a*x+b

class Simulation(object):
    def __init__(self):
        pygame.init()

        self.screen = pygame.display.set_mode([SCREEN_SIZE, SCREEN_SIZE])
        self.clock = pygame.time.Clock()

        self.coordinates = []
        self.velocities = []
        self.wall_collides = []
        self.particle_collides = []
        self.collision_times = []
        self.collision_form = []

        self.init_disks()
        self.init_wall_collides()
        self.init_particle_collides()

        if WRITE_PDF:
            self.velocity_dist = open("velo.dat", "w")
            self.pp = PdfPages("velo_dist.pdf")
        elif PRESSURE:
            self.temp_points = []
            self.pressure_points = []
            self.pressure_dat = open("pressure.dat", "w")

    def reset(self):
        self.coordinates = []
        self.velocities = []
        self.wall_collides = []
        self.particle_collides = []
        self.collision_times = []
        self.collision_form = []

        self.init_disks()
        self.init_wall_collides()
        self.init_particle_collides()

    def run(self):
        if PRESSURE:
            for i in range(3):                
                print(i)
                self.do_simulation(PRESSURE_LENGTH)
                global MAX_VELOCITY
                MAX_VELOCITY += 1
                self.reset()
            self.write_pressure_graph()
            self.pressure_dat.close()
        else:
            self.do_simulation()
    
    def do_simulation(self, length = INFINITY):
        collide = True
        t = 0
        v_dist_interval = 0
        pressure_dat_interval = 0
        pressure_per_interval = 10
        pressure_dat = []

        try:
            while t < length:
                if t >= v_dist_interval and WRITE_PDF:
                    self.write_velocities()
                    v_dist_interval += 50
                elif t >= pressure_dat_interval and PRESSURE:
                    pressure_dat.append(pressure_per_interval)
                    pressure_dat_interval += 10

                #self.print_velocities()
                #print(t)
                if collide:
                    self.find_collision_times()
                    next_collide = min(self.collision_times)
                    particle = self.collision_times.index(next_collide)
                    collide = False
                    # if next_collide < t:
                    # print("nc: "+str(next_collide))
                    # print("t: "+str(t))
                    # print(self.velocities)
                    # print(self.coordinates)
                    # print(self.particle_collides)
                    # print(self.wall_collides)
                    # print(self.collision_times)
                    # print(self.collision_form)

                self.clock.tick(CPU_TICKS)
                self.screen.fill(WHITE)

                for i in self.coordinates:
                    pygame.draw.circle(self.screen, BLUE, [int(round(j)) for j in i], RADIUS)

                pygame.display.flip()

                if abs(next_collide - t) <= DT:
                    self.do_tick(abs(next_collide - t))

                    if self.collision_form[particle][0] == WALL_COLLISION:
                        direction = self.collision_form[particle][1]
                        vx = self.velocities[particle][0]
                        vy = self.velocities[particle][1]

                        if direction == X_DIRECTION:
                            self.velocities[particle] = [-vx, vy]
                            if PRESSURE:
                                pressure_per_interval += abs(vx)/SCREEN_SIZE
                        elif direction == Y_DIRECTION:
                            self.velocities[particle] = [vx, -vy]
                            if PRESSURE:
                                pressure_per_interval += abs(vy)/SCREEN_SIZE

                        self.calc_new_wall_collides(particle, t, direction)
                        self.calc_new_particle_collides(t, particle)
                    else:
                        other_particle = self.collision_form[particle][0]

                        r12 = vec_minus(self.coordinates[particle], self.coordinates[other_particle])

                        v1 = self.velocities[particle]
                        v1p = rescale(r12, dot(v1, r12)/dot(r12, r12))
                        v1s = vec_minus(v1, v1p)

                        v2 = self.velocities[other_particle]
                        v2p = rescale(r12, dot(v2, r12)/dot(r12, r12))
                        v2s = vec_minus(v2, v2p)

                        v1pnew = v2p
                        v2pnew = v1p
                        v1new = vec_plus(v1pnew, v1s)
                        v2new = vec_plus(v2pnew, v2s)

                        self.velocities[particle] = v1new
                        self.velocities[other_particle] = v2new

                        self.calc_new_wall_collides(particle, t)
                        self.calc_new_wall_collides(other_particle, t)
                        self.calc_new_particle_collides(t, particle, other_particle)

                    collide = True
                    t += abs(next_collide - t)
                else:
                    self.do_tick()
                    t += DT
        except KeyboardInterrupt:
            pass
        if WRITE_PDF:
            self.pp.close()
            self.velocity_dist.close()
        elif PRESSURE:
            T = sum(dot(i, i) for i in self.velocities)/N
            pressure = sum(i for i in pressure_dat)/len(pressure_dat)
            self.temp_points.append(T)
            self.pressure_points.append(pressure)
            self.pressure_dat.write(str(T)+" "+str(pressure)+"\n")

    def write_velocities(self):
        # abs_velocity = []
        # data = []

        # for i in self.velocities:
        #     abs_velocity.append(dot(i, i))

        # for i in range(0, 60, 5):
        #     counter = 0
        #     abs_velocity_copy = []
        #     for j in abs_velocity:
        #         if j <= i/6:
        #             counter += 1
        #         else:
        #             abs_velocity_copy.append(j)
        #     abs_velocity = abs_velocity_copy[:]
        #     data.append(counter)

        for i in self.velocities:
            self.velocity_dist.write(str(dot(i, i))+ " ")
        self.velocity_dist.write("\n")

        P.xlabel("velocity v")
        P.ylabel("Propability P")
        P.title("Velocity distribution")
        n, bins, patches = P.hist(P.array([dot(i, i) for i in self.velocities]), [i/3. for i in range(30)], normed=1, histtype="bar")
        self.pp.savefig()
        P.figure()

    def write_pressure_graph(self):
        P.xlabel("Temperature")
        P.ylabel("Pressure")
        P.title("Pressure per Temperature")
        P.axis([0.0,self.temp_points[-1]+10., 0.0,self.pressure_points[-1]+1])
        ax = P.gca()
        ax.set_autoscale_on(False)

        popt,pcov = curve_fit(fit,self.temp_points,self.pressure_points)
        y_fit = [popt[0]*x+popt[1] for x in self.temp_points]
        y_fit.insert(0, popt[1])
        fit_x_points = self.temp_points[:]
        fit_x_points.insert(0,0.)

        pp = PdfPages(str(SCREEN_SIZE)+"_"+str(N)+".pdf")
        P.plot(self.temp_points, self.pressure_points, "o", fit_x_points, y_fit, "--")
        #P.savefig()
        #P.plot(
        #pp.savefig()
        #print(self.pressure_points)
        #print(y_fit)
        #print(fit_x_points)
        #print(y_fit)
        pp.savefig()
        pp.close()

    def do_tick(self, dt = DT):
        for i in range(len(self.coordinates)):
            self.coordinates[i][0] += dt * self.velocities[i][0]
            self.coordinates[i][1] += dt * self.velocities[i][1]

    def init_disks(self):
        i = 0
        while i < N:
            x = random.randint(RADIUS + 1, SCREEN_SIZE - (RADIUS + 1))
            y = random.randint(RADIUS + 1, SCREEN_SIZE - (RADIUS + 1))

            try:
                for j in self.coordinates:
                    if norm(vec_minus([x, y], j)) <= 2*RADIUS:
                        raise Exception
            except:
                continue

            if EQUILIBRIUM:
                a = random.random()
                b = random.random()
                vx = MAX_VELOCITY*math.sqrt(-2*math.log(1-a))*math.cos(2*math.pi*b)
                vy = MAX_VELOCITY*math.sqrt(-2*math.log(1-a))*math.sin(2*math.pi*b)
            else:
                direction = random.randint(0, 3)
                if direction == 0:
                    vx = MAX_VELOCITY
                    vy = 0
                elif direction == 1:
                    vx = 0
                    vy = MAX_VELOCITY
                elif direction == 2:
                    vx = -MAX_VELOCITY
                    vy = 0
                elif direction == 3:
                    vx = 0
                    vy = -MAX_VELOCITY
                    # vx = random.randint(-MAX_VELOCITY, MAX_VELOCITY)
                    # vy = random.randint(-MAX_VELOCITY, MAX_VELOCITY)

            self.coordinates.append([x, y])
            self.velocities.append([vx, vy])
            i += 1
                
    def init_wall_collides(self):
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

    def calc_new_wall_collides(self, particle, t, direction = None):
        if direction == 0 or (not direction):
            if direction:
                ty = self.wall_collides[particle][1]

            x = self.coordinates[particle][0]
            vx = self.velocities[particle][0]

            if vx > 0:
                tx = (SCREEN_SIZE - (RADIUS + x)) / vx
                tx += t
            elif vx < 0:
                tx = (RADIUS - x) / vx
                tx += t
            else:
                tx = INFINITY
        if direction == 1 or (not direction):
            if direction:
                tx = self.wall_collides[particle][0]

            y = self.coordinates[particle][1]
            vy = self.velocities[particle][1]

            if vy > 0:
                ty = (SCREEN_SIZE - (RADIUS + y)) / vy
                ty += t
            elif vy < 0:
                ty = (RADIUS - y) / vy
                ty += t
            else:
                ty = INFINITY

        self.wall_collides[particle] = [tx, ty]

    def init_particle_collides(self):
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

    def calc_new_particle_collides(self, t, particle1, particle2 = None):
        times1 = []
        times2 = []
        calc1 = True
        calc2 = True

        r1 = self.coordinates[particle1]
        v1 = self.velocities[particle1]
        if particle2:
            r2 = self.coordinates[particle2]
            v2 = self.velocities[particle2]

        #print("particle1 "+str(particle1))
        #print("particle2 "+str(particle2))
        #print(self.particle_collides)
        for j in range(len(self.coordinates)):
            if j == particle1:
                times1.append(INFINITY)
                calc1 = False
            elif particle2 is not None and j == particle2:
                times2.append(INFINITY)
                calc2 = False

            rj = self.coordinates[j]
            vj = self.velocities[j]

            if calc1:
                r1j = vec_minus(r1, rj)
                v1j = vec_minus(v1, vj)
                d = (dot(r1j, v1j)**2 - dot(v1j, v1j) * (dot(r1j, r1j) - 4 * RADIUS**2))
                if d > 0:
                    tcoll = (-dot(r1j, v1j)-math.sqrt(d))/dot(v1j, v1j)
                    if tcoll > 0:
                        times1.append(t + tcoll)
                    else:
                        times1.append(INFINITY)
                else:
                    times1.append(INFINITY)
            if particle2 and calc2:
                r2j = vec_minus(r2, rj)
                v2j = vec_minus(v2, vj)
                d = (dot(r2j, v2j)**2 - dot(v2j, v2j) * (dot(r2j, r2j) - 4 * RADIUS**2))
                if d > 0:
                    tcoll = (-dot(r2j, v2j)-math.sqrt(d))/dot(v2j, v2j)
                    if tcoll > 0:
                        times2.append(t + tcoll)
                    else:
                        times2.append(INFINITY)
                else:
                    times2.append(INFINITY)
            calc1 = True
            calc2 = True

        #print("times1: "+str(times1))
        #print("times2: "+str(times2))
    
        self.particle_collides[particle1] = times1
        if particle2: 
            self.particle_collides[particle2] = times2

        for i in range(len(self.coordinates)):
            self.particle_collides[i][particle1] = self.particle_collides[particle1][i]
            if particle2: 
                self.particle_collides[i][particle2] = self.particle_collides[particle2][i]
        #print(self.particle_collides)

    def find_collision_times(self):
        self.collision_times = []
        self.collision_form = []
        for i in range(len(self.coordinates)):
            t_wall = self.wall_collides[i]
            if t_wall[0] < t_wall[1]:
                t_wall = t_wall[0]
                direction = X_DIRECTION
            else:
                t_wall = t_wall[1]
                direction = Y_DIRECTION

            t_particle = min(self.particle_collides[i])

            #print(t_wall)
            #print(t_particle)
            if t_wall <= t_particle:
                #print(1)
                self.collision_times.append(t_wall)
                self.collision_form.append([WALL_COLLISION, direction])
            else:
                #print(2)
                self.collision_times.append(t_particle)
                self.collision_form.append([self.particle_collides[i].index(t_particle)])


if __name__ == "__main__":
    sim = Simulation()
    sim.run()
            
