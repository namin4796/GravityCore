import sys
import os
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# load c++ engine
sys.path.append(os.path.join(os.path.dirname(__file__), '../build'))
try:
    import gravity_core
    print("SUCCESS: Loaded C++ gravity_core engine!")
except ImportError:
    print("ERROR: Could not load 'gravity_core'. Did you run 'cmake --build .' inside build/?")
    sys.exit(1)

# -- config --
N = 1000 #Number of stars
STEPS_PER_FRAME = 20 # how many physics steps to take

# initialize universe
sim = gravity_core.Universe(N)

# create a galaxy using polar co-ordinates
angle = [random.uniform(0, 6.28) for _ in range(N)]
radius = [random.uniform(10, 100) for _ in range(N)]

x = [r * np.cos(a) for r, a in zip(radius, angle)]
y = [r * np.sin(a) for r, a in zip(radius, angle)]
mass = [random.uniform(1.0, 5.0) for _ in range(N)]

# add a SMB in the centre
x[0], y[0] = 0, 0
mass[0] = 50000.0

# initial velocity
vx = []
vy = []
for i in range(N):
    r = radius[i]
    if r < 1:
        vx.append(0)
        vy.append(0)
        continue

    v_orbit = np.sqrt(1.0 * mass[0] / r) * 0.7

    # perpendicular vel. to radius
    vx.append(-v_orbit * np.sin(angle[i]))
    vy.append(v_orbit * np.cos(angle[i]))

# set positions and mass
sim.set_state(x, y, vx, vy, mass)

# -- visualization loop --
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_facecolor('black')
ax.set_xlim(-150, 150)
ax.set_ylim(-150, 150)

# scatter plot for stars
particles = ax.scatter([], [], c='white', s=2, alpha=0.8)
# the SMB
center = ax.scatter([], [], c='red', s=5)

def update(frame):
    # run c++ physics
    for _ in range(STEPS_PER_FRAME):
        sim.step()

    # get data from c++; instant for 1000 particles
    # convert C++ vector to a python list/array for matplotlib
    new_x = np.array(sim.get_x())
    new_y = np.array(sim.get_y())

    # print first particle's position
    # if it prints same number every time, the c++ is frozen.
    if frame % 10 == 0:
        print(f"Frame {frame}: Particle 1 is at ({new_x[1]:.2f}, {new_y[1]:.2f})")

    # update plot
    data = np.c_[new_x, new_y]
    particles.set_offsets(data)
    center.set_offsets(np.c_[[new_x[0]], [new_y[0]]])
    return particles, center

print("Starting Animation... Close window to exit.")
ani = FuncAnimation(fig, update, frames=200, interval=20, blit=False)
plt.show()
