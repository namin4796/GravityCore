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
    show_tree = False

    print("Controls: Press 'T' to toggle Quadtree visualization. Press 'ESC' to quit.")


    while running:
        # event handler
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_ESCAPE:
                    running = False
                if event.key == pygame.K_t:
                    show_tree = not show_tree

        engine.step()

        # opengl rendering
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        # draw quadtree
        if show_tree:
            rects = engine.get_tree_rects()
            if len(rects) > 0:
                # convert the flat list into struct. Nx3 array
                rects_np = np.array(rects, dtype=np.float64).reshape(-1, 3)

                x = rects_np[:, 0]
                y = rects_np[:, 1]
                s = rects_np[:, 2]

                #vectorized generation of line
                lines = np.empty((len(rects_np), 8, 2), dtype=np.float64)

                lines[:, 0, 0] = x;     lines[:, 0, 1] = y  # bottom-left to 
                lines[:, 1, 0] = x + s; lines[:, 1, 1] = y  # bottom-right
                
                lines[:, 2, 0] = x + s; lines[:, 2, 1] = y  # bottom-right to
                lines[:, 3, 0] = x + s; lines[:, 3, 1] = y + s # top-right

                lines[:, 4, 0] = x + s; lines[:, 4, 1] = y + s   # Top-right to
                lines[:, 5, 0] = x;     lines[:, 5, 1] = y + s   # Top-left
                
                lines[:, 6, 0] = x;     lines[:, 6, 1] = y + s   # Top-left to
                lines[:, 7, 0] = x;     lines[:, 7, 1] = y       # Bottom-left

                lines_gl = np.ascontiguousarray(lines.reshape(-1, 2))

                # draw all line at once
                glEnableClientState(GL_VERTEX_ARRAY)
                glVertexPointer(2, GL_DOUBLE, 0, lines_gl)
                glColor3f(0.15, 0.15, 0.15)
                glDrawArrays(GL_LINES, 0, len(lines_gl))
                glDisableClientState(GL_VERTEX_ARRAY)
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
