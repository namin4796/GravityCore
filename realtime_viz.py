import sys
import os
import json
import pygame
from pygame.locals import *
from OpenGL.GL import *
from OpenGL.GLU import *
import numpy as np

script_dir = os.path.dirname(os.path.abspath(__file__))
build_path = os.path.join(script_dir, "build")
sys.path.insert(0, build_path)

try:
    import gravity_core
except ImportError as e:
    print(f"CRITICAL ERROR: could not find 'gravity_core' in {build_path}")
    print(f"Original error: {e}")
    sys.exit(1)


def load_config(filename):
    try:
        with open(filename, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        print(f"CRITICAL ERROR: Could not find config file '{filename}'")
        sys.exit(1)

def main():
    if len(sys.argv) < 2:
        print("Usage: python realtime_viz.py <config.json>")
        sys.exit(1)

    config_file = sys.argv[1]
    print(f"Loading parameters from {config_file}...")
    config = load_config(config_file)

    num_stars = config.get("NUM_STARS", 5000)
    dt = config.get("DT", 0.01)
    r_scale = config.get("R_SCALE", 100.0)
    rho_0 = config.get("RHO_0", 0.01)
    G = config.get("G", 2.0)

    #galaxy specs
    gal_radius = config.get("GAL_RAD", 800.0)
    central_mass = config.get("MASS_BH", 100000.0)

    #initialize engine
    pygame.init()
    WIDTH, HEIGHT = 1200, 1000
    pygame.display.set_caption(f"GravityCore: {config_file}")
    pygame.display.set_mode((WIDTH, HEIGHT), DOUBLEBUF | OPENGL)

    #setup opengl camera
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(45, (WIDTH / HEIGHT), 0.1, 50000.0)
    glMatrixMode(GL_MODELVIEW)

    print(f"Initializing Physics with {num_stars} bodies...")
    engine = gravity_core.Universe(num_stars, r_scale, rho_0, G, dt)

    #initialize with a disk galaxy setup
    engine.init_galaxy(gal_radius, central_mass)

    # generating realistic colors
    # core stars brighter (yellow/white), outer stars blue
    pos = engine.get_positions()
    radii = np.sqrt(pos[:, 0]**2 + pos[:, 1]**2)
    max_radius = np.max(radii) if np.max(radii) > 0 else 1.0
    norm_radii = np.clip(radii / (max_radius * 0.5), 0, 1)

    colors = np.zeros((num_stars, 3), dtype=np.float32)
    colors[:, 0] = 1.0 - norm_radii * 0.6
    colors[:, 1] = 0.9 - norm_radii * 0.4
    colors[:, 2] = 0.7 + norm_radii * 0.3
    colors_gl = np.ascontiguousarray(colors)

    # camera state
    cam_x, cam_y, cam_z = 0.0, 0.0, -2500.0
    dragging = False
    last_mouse = (0, 0)
    show_tree = False

    clock = pygame.time.Clock()
    running = True

    print("Controls: Press 'T' to toggle Quadtree visualization. Press 'ESC' to quit.")


    while running:
        # event handler
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
            elif event.type == pygame.KEYDOWN:
                if event.key == pygame.K_ESCAPE:
                    running = False
                if event.key == pygame.K_t:
                    show_tree = not show_tree
            
            #mouse control
            elif event.type == pygame.MOUSEBUTTONDOWN:
                if event.button == 4:
                    cam_z += 200
                elif event.button == 5:
                    cam_z -= 200
                elif event.button == 1:
                    dragging = True
                    last_mouse = event.pos
            elif event.type == pygame.MOUSEBUTTONUP:
                if event.button == 1:
                    dragging = False
            elif event.type == pygame.MOUSEMOTION:
                if dragging:
                    dx, dy = event.pos[0] - last_mouse[0], event.pos[1] - last_mouse[1]
                    scale = abs(cam_z) / 1000.0
                    cam_x += dx * scale
                    cam_y -= dy * scale
                    last_mouse = event.pos

        engine.step()

        # opengl rendering
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glLoadIdentity()

        glTranslatef(cam_x, cam_y, cam_z) # zoom out to see the disc
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
                glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
                glEnable( GL_BLEND )
                glColor4f(0.5, 1.0, 0.5, 0.25)
                glDrawArrays(GL_LINES, 0, len(lines_gl))
                glDisableClientState(GL_VERTEX_ARRAY)
        # fetch data from C++
        pos = engine.get_positions()

        pos_gl = np.ascontiguousarray(pos, dtype=np.float64)

        # optimize OpenGL rendering
        glEnableClientState(GL_VERTEX_ARRAY)
        glEnableClientState(GL_COLOR_ARRAY)
        # using double since C++ uses double
        glVertexPointer(2, GL_DOUBLE, 0, pos_gl)
        glColorPointer(3, GL_FLOAT, 0, colors_gl)

        glPointSize(2.0)

        current_num_stars = len(pos_gl) # dynamically checking array size
        glDrawArrays(GL_POINTS, 0, current_num_stars)

        glDisableClientState(GL_COLOR_ARRAY)
        glDisableClientState(GL_VERTEX_ARRAY)

        pygame.display.flip()
        clock.tick(60) # 60 FPS

        print(f"FPS: {clock.get_fps():.2f} | Z: {cam_z:.0f}", end='\r')
    
    pygame.quit()

if __name__ == "__main__":
    main()
