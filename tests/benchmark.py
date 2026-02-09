import sys
import os
import time
import argparse
import numpy as np

# load engine
sys.path.append(os.path.join(os.path.dirname(__file__), '../build'))
try:
    import gravity_core
except ImportError:
    print("Error: Could not load gravity_core. Compile first!")
    sys.exit(-1)

def run_benchmark(n_stars, steps):
    print(f"--BENCHMARKING N={n_stars}---")

    # use uge radius to ensure deep tree structure
    sim = gravity_core.Universe(n_stars, 100.0, 0.001, 4.3e-6, 0.01)

    # random data
    x = np.random.uniform(-200, 200, n_stars)
    y = np.random.uniform(-200, 200, n_stars)
    mass = np.random.uniform(1.0, 10.0, n_stars)
    vx = np.zeros(n_stars)
    vy = np.zeros(n_stars)

    sim.set_state(x, y, vx, vy, mass)

    # warmup (let the OS allocate memory pages)
    print("Warming up ...")
    for _ in range(10):
        sim.step()

    # race
    print(f"Running {steps} steps...", end="", flush=True)
    start_time = time.time()

    for _ in range(steps):
        sim.step()

    end_time = time.time()
    duration = end_time - start_time

    # results
    sps = steps / duration
    interactions = n_stars * steps
    # O(N log N) estimation is hard to quantify exactly, so we track raw interactions

    print(f"DONE!")
    print(f"Time: {duration:.4f} seconds")
    print(f"Speed: {sps:.2f} steps/second")
    print(f"-----------------------------")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, default=5000)
    parser.add_argument("--STEPS", type=int, default=100)
    args = parser.parse_args()

    run_benchmark(args.N, args.STEPS)
