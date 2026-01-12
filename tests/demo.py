import sys
import os
import time
import random
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), '../build'))

import gravity_core

N = 2000
print(f"Initializing Universe with {N} particles...")

#create engine
sim = gravity_core.Universe(N)

#setup data (star cluster) random positions
x = [random.uniform(-100, 100) for _ in range(N)]
y = [random.uniform(-100, 100) for _ in range(N)]
mass = [random.uniform(1.0, 10.0) for _ in range(N)]

#make first particle a supermassive black hole
x[0], y[0] = 0, 0
mass[0] = 1000000.0

#load data into c++
sim.set_state(x, y, mass)

print("Starting Simulation (C++)...")
start = time.time()

steps = 100
for i in range(steps):
    sim.step()

duration = time.time() - start
print(f"Done! {steps} for {N} particles took {duration: .4f}s")
print(f"Performance: {(steps*N*N)/duration/1e6:.2f} Million interactions/sec")

#python comparison
# a pure python double-loop for 2000 particles takes ~0.5s PER STEP.
# 100 steps would take ~50 seconds.
print(f"Estimated pure Python time: {steps * 0.5:.2f}s")
print(f"Speedup Factor: ~{ (steps * 0.5) /duration : .1f}x")
