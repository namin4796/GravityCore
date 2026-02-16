import sys
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
build_path = os.path.join(script_dir, "build")
sys.path.append(build_path)

import pygame
from pygame.locals import *
from OpenGL.GL import *
from OpenGL.GLU import *
import numpy as np

try:
    import gravity_core
except ImportError as e:
    print(f"CRITICAL ERROR: could not find 'gravity_core' in {build_path}")
    print(f"Original error: {e}")
    sys.exit(1)

WIDTH, HEIGHT = 1200, 1000
NUM_STARS = 5000

def main():
    pygame.init()
    pygame.display.set_caption("")
    pygame.display.set_mode((WIDTH, HEIGHT), DOUBLEBUF | OPENGL)

    #setup opengl camera
    gluPerspective(45, (WIDTH / HEIGHT), 0.1, 10000.0)
    glTranslatef(0.0, 0.0, -2000) # zoom out to see the disc

    print(f"Initializing Physics with {NUM_STARS} bodies...")
    engine = gravity_core.Universe(NUM_STARS, 100.0, 0.01, 5.0, 0.1)

    #initialize with a disk galaxy setup
    engine.init_galaxy(800.0, 100000.0)

    clock = pygame.time.Clock()
    running = True

    while running:
        # event handler
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_ESCAPE:
                    running = False

        engine.step()

        # opengl rendering
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        # fetch data from C++
        pos = engine.get_positions()

        pos_gl = np.ascontiguousarray(pos, dtype=np.float64)

        # optimize OpenGL rendering
        glEnableClientState(GL_VERTEX_ARRAY)
        # using double since C++ uses double
        glVertexPointer(2, GL_DOUBLE, 0, pos_gl)

        glPointSize(2.0)

        glColor3f(0.4, 0.7, 1.0)
        glDrawArrays(GL_POINTS, 0, NUM_STARS)

        glDisableClientState(GL_VERTEX_ARRAY)

        pygame.display.flip()
        clock.tick(60) # 60 FPS

        print(f"FPS: {clock.get_fps():.2f}", end='\r')
    
    pygame.quit()

if __name__ == "__main__":
    main()
