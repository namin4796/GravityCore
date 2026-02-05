import argparse
import sys
import os
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

parser = argparse.ArgumentParser(description="Run N-body Simulation")
parser.add_argument("--N_STARS", type=int, default=1000, help="Number of Stars")
parser.add_argument("--R_SCALE", type=float, default=100.0, help="NFW Scale Radius")
parser.add_argument("--RHO_0", type=float, default=0.001, help="Dark Matter Density")
parser.add_argument("--DT", type=float, default=0.01, help="Time step per physics tick")
parser.add_argument("--STEPS", type=int, default=1, help="Physics steps per render frame")
parser.add_argument("--DEBUG_TREE", action="store_true", help="Visualize Barnes0-Hut Quadtree")
args = parser.parse_args()

def get_total_mass(r, m_bh, R_SCALE, RHO_0):
    x = r / R_SCALE
    m_NFW = 4.0 * np.pi * RHO_0 * (R_SCALE**3)
    m_Halo = m_NFW * (np.log(1+x)-(x/(1+x)))

    return (m_bh + m_Halo)

# load c++ engine
sys.path.append(os.path.join(os.path.dirname(__file__), '../build'))
try:
    import gravity_core
    print("SUCCESS: Loaded C++ gravity_core engine!")
except ImportError:
    print("ERROR: Could not load 'gravity_core'. Did you run 'cmake --build .' inside build/?")
    sys.exit(1)

# astrophysical constants
# unit systems: [kpc, solar mass, km/s]
G_GALACTIC = 4.302e-6
DT_SAFE = 0.000001
# milky way parameters
# scale radius ~ 15 - 20 kpc
# R_SCALE_MW = 20.0
# dark matter density ~ 0.008 M_sun/pc^3 -> convert M_sun/kpc^3
# 0.008 * 10^9 = 8,000,000 M_sun / kpc^3
# RHO_0_MW = 8.0e6
MASS_BH = 4.0e6 # 4 Million Solar masses

# -- config --
N = args.N_STARS  #Number of stars
STEPS_PER_FRAME = args.STEPS # how many physics steps to take

# initialize universe
sim = gravity_core.Universe(N, args.R_SCALE, args.RHO_0, G_GALACTIC, args.DT)

# create a galaxy using polar co-ordinates
angle = [random.uniform(0, 2*np.pi) for _ in range(N)]
radius = [random.uniform(1.0, 50.0) for _ in range(N)]

x = [r * np.cos(a) for r, a in zip(radius, angle)]
y = [r * np.sin(a) for r, a in zip(radius, angle)]
mass = [random.uniform(10.0, 50.0) for _ in range(N)]

# add a SMB in the centre
x[0], y[0] = 0, 0
mass[0] = MASS_BH

# initial velocity
vx = []
vy = []
for i in range(N):
    r = radius[i]
    if r < 1:
        vx.append(0)
        vy.append(0)
        continue

    m_Total = get_total_mass(r, mass[0], args.R_SCALE, args.RHO_0)
    v_orbit = np.sqrt(G_GALACTIC * m_Total / r)

    # small noise to v
    v_orbit *= random.uniform(0.95, 1.05) #slight variation

    # perpendicular vel. to radius (tangential vel.)
    vx.append(-v_orbit * np.sin(angle[i]))
    vy.append(v_orbit * np.cos(angle[i]))

# set positions and mass
sim.set_state(x, y, vx, vy, mass)

# -- visualization loop --
fig, ax = plt.subplots(figsize=(10, 10))
ax.set_facecolor('black')
ax.set_xlim(-80, 80)
ax.set_ylim(-80, 80)

# scatter plot for stars
particles = ax.scatter([], [], c='cyan', s=1, alpha=0.8)
# the SMB
center = ax.scatter([], [], c='yellow', s=15, zorder=10)

# collection for tree boxes
rect_collection = PatchCollection([], edgecolors='green', facecolors='none', alpha=0.2, linewidths=0.5)
ax.add_collection(rect_collection)
                                 

def update(frame):
    # run c++ physics
    for _ in range(STEPS_PER_FRAME):
        sim.step()

    # get data from c++; instant for 1000 particles
    # convert C++ vector to a python list/array for matplotlib
    new_x = np.array(sim.get_x())
    new_y = np.array(sim.get_y())


    # update plot
    data = np.c_[new_x, new_y]
    particles.set_offsets(data)
    center.set_offsets(np.c_[[new_x[0]], [new_y[0]]])
    
    # DEBUG: Draw Qaudtree
    if args.DEBUG_TREE:
        rect_data = sim.get_tree_rects()
        patches = []
        # limit to first 500 boxes to prevent lag if tree is huge
        limit = min(len(rect_data), 1500)
        for i in range(0, limit, 3):
            rx, ry, size = rect_data[i], rect_data[i+1], rect_data[i+2]
            patches.append(Rectangle((rx, ry), size, size))

        rect_collection.set_paths(patches)

    return particles, center, rect_collection

print("Starting Animation... Close window to exit.")
ani = FuncAnimation(fig, update, frames=200, interval=1, blit=False)
plt.title(f"Galaxy Simulation (N={N} - Quadtree Viz: {'ON' if args.DEBUG_TREE else 'OFF'}")
plt.xlabel("R / kpc")
plt.ylabel("R / kpc")
plt.show()
