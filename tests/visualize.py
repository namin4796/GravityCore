import argparse
import sys
import os
import json
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

parser = argparse.ArgumentParser(description="Run N-body Simulation")
parser.add_argument("--config", type=str, default="default_config.json", help="Path to JSON config file.")
parser.add_argument("--DEBUG_TREE", action="store_true", help="Visualize Barnes0-Hut Quadtree")
parser.add_argument("--N_STARS", type=int, default=1000, help="Number of Stars")
parser.add_argument("--DT", type=float, default=0.01, help="Time step per physics tick")
args = parser.parse_args()

# Load Defaults from Config
config = {}
if os.path.exists(args.config):
    with open(args.config, 'r') as f:
        config = json.load(f)
else:
    print(f"WARNING: Config file {args.config} not found. Using hardcoded defaults.")

# parsing overriding parameters
N = args.N_STARS if args.N_STARS else config.get("N_STARS", 1000)
DT = args.DT if args.DT else config.get("DT", 0.01)
R_SCALE = config.get("R_SCALE", 100.0)
RHO_0 = config.get("RHO_0", 0.001)
STEPS_PER_FRAME = config.get("STEPS", 1)
THETA = config.get("THETA", 1.0) # vel. multiplier
MASS_BH = config.get("MASS_BH", 4.0e6)

# load c++ engine
sys.path.append(os.path.join(os.path.dirname(__file__), '../build'))
try:
    import gravity_core
    print("SUCCESS: Loaded C++ gravity_core engine!")
except ImportError:
    print("ERROR: Could not load 'gravity_core'. Did you run 'cmake --build .' inside build/?")
    sys.exit(1)

# helper function to get total mass
def get_total_mass(r, m_bh):
    x = r / R_SCALE
    m_NFW = 4.0 * np.pi * RHO_0 * (R_SCALE**3)
    m_Halo = m_NFW * (np.log(1+x)-(x/(1+x)))

    return (m_bh + m_Halo)
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
#MASS_BH = 4.0e6 # 4 Million Solar masses


# initialize universe
sim = gravity_core.Universe(N, R_SCALE, RHO_0, G_GALACTIC, DT)

# create a galaxy using polar co-ordinates
angle = [random.uniform(0, 2*np.pi) for _ in range(N)]
radius = [random.uniform(2.0, 60.0) for _ in range(N)]

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
    if i == 0:
        vx.append(0)
        vy.append(0)
        continue

    m_Total = get_total_mass(r, mass[0])
    v_orbit = np.sqrt(G_GALACTIC * m_Total / r)

    # small noise to v
    v_orbit *= THETA #slight variation

    # adding randomness to direction
    v_radial = 0.0
    if THETA < 0.5:
        v_radial = v_orbit*0.5

    # perpendicular vel. to radius (tangential vel.)
    vx.append(-v_orbit * np.sin(angle[i]) + v_radial * np.cos(angle[i]))
    vy.append(v_orbit * np.cos(angle[i])+ v_radial * np.sin(angle[i]))

# set positions and mass
sim.set_state(x, y, vx, vy, mass)

# -- visualization loop --
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_facecolor('black')
ax.set_xlim(-80, 80)
ax.set_ylim(-80, 80)

# Init color map
cmap = plt.get_cmap('inferno')

# scatter plot for stars
particles = ax.scatter([], [], c=[], cmap=cmap, s=2, alpha=0.8)
# the SMB
center = ax.scatter([], [], c='white', s=20, zorder=10)

# collection for tree boxes
rect_collection = PatchCollection([], edgecolors='green', facecolors='none', alpha=0.15, linewidths=0.5)
ax.add_collection(rect_collection)
                                 

def update(frame):
    # run c++ physics
    for _ in range(STEPS_PER_FRAME):
        sim.step()

    # get data from c++; instant for 1000 particles
    # convert C++ vector to a python list/array for matplotlib
    new_x = np.array(sim.get_x())
    new_y = np.array(sim.get_y())
    new_vx = np.array(sim.get_vx())
    new_vy = np.array(sim.get_vy())


    # update plot
    data = np.c_[new_x, new_y]
    particles.set_offsets(data)
    center.set_offsets(np.c_[[new_x[0]], [new_y[0]]])
    
    # calculate magnitude of velocity
    v_mag = np.sqrt(new_vx**2 + new_vy**2)

    # Normalize for color map:
    # We clip the visualization at v=200 km/s (or whatever unit)
    # to prevent one crazy star from making everyone else black.
    # Adjust this V_MAX if stars are too dark or too washed out.
    V_MAX = 0.8

    # Skip index 0 (Black Hole) for coloring to avoid skewing
    colors = v_mag
    colors[0] = 0 # Hide BH temp

    particles.set_array(colors)
    particles.set_clim(0, V_MAX) # Set the range for the colormap

    # DEBUG: Draw Qaudtree
    if args.DEBUG_TREE:
        rect_data = sim.get_tree_rects()
        patches = []
        # limit to first 500 boxes to prevent lag if tree is huge
        limit = min(len(rect_data), 5000)
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
