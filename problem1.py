import pygame
import random
import math

WHITE = (255, 255, 255)
BLUE =  (  0,   0, 255)
GREEN = (  0, 255,   0)
RED =   (255,   0,   0)

N = 100
RADIUS = 8
MAX_VELOCITY = 5
SCREEN_SIZE = 400
CPU_TICKS = 20

class Simulation(object):
    def __init__(self):
        pygame.init()

        self.screen = pygame.display.set_mode([SCREEN_SIZE, SCREEN_SIZE])
        self.clock = pygame.time.Clock()

        self.coordinates = []
        self.velocities = []

        self.init_disks()

    def run(self):

        while True:
            self.clock.tick(CPU_TICKS)
            self.screen.fill(WHITE)

            for i in self.coordinates:
                pygame.draw.circle(self.screen, BLUE, i, RADIUS)

            pygame.display.flip()
            self.do_tick()

    def do_tick(self):
        dt = 1

        for i in range(len(self.coordinates)):
            self.coordinates[i][0] += dt * self.velocities[i][0]
            self.coordinates[i][1] += dt * self.velocities[i][1]

    def init_disks(self):
        i = 0
        while i < N:
            x = random.randint(RADIUS + 1, SCREEN_SIZE - (RADIUS + 1))
            y = random.randint(RADIUS + 1, SCREEN_SIZE - (RADIUS + 1))

            for j in self.coordinates:
                if math.sqrt((x - j[0])**2 + (y - j[1])**2) <= RADIUS:
                    continue

            vx = random.randint(-MAX_VELOCITY, MAX_VELOCITY)
            vy = random.randint(-MAX_VELOCITY, MAX_VELOCITY)

            self.coordinates.append([x, y])
            self.velocities.append([vx, vy])

            i += 1 
                

if __name__ == "__main__":
    sim = Simulation()
    sim.run()
            
