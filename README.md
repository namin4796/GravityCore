# GravityCore: High-Performance N-Body Physics Engine

![Language](https://img.shields.io/badge/language-C%2B%2B17%20%7C%20Python-blue)
![Build](https://img.shields.io/badge/build-CMake-green)
![Tests](https://img.shields.io/badge/tests-Catch2-success)
![License](https://img.shields.io/badge/license-MIT-lightgrey)

**GravityCore** is a hybrid physics simulation engine capable of calculating gravitational interactions for thousands of bodies in real-time at 60 FPS. It leverages **Modern C++ (C++17)** for high-performance computing and **Hardware-Accelerated OpenGL** for visualization via a Zero-Copy Python memory bridge.

## Key Architecture

### 1. $O(N \log N)$ Barnes-Hut Spatial Partitioning
Instead of computationally expensive $O(N^2)$ all-to-all force calculations, the engine implements a **Quadtree**. Distant particle clusters are approximated as single centers of mass, drastically reducing algorithmic complexity and allowing for massive scale.

### 2. Zero-Copy Memory Bridge
Data is never copied between C++ and Python during the simulation loop. The C++ Structure of Arrays (SoA) memory is exposed directly to Python as contiguous `numpy` arrays via `pybind11`, which are fed instantly to the GPU using OpenGL Vertex Arrays (`glDrawArrays`).

### 3. Astrophysical Realism
The engine goes beyond simple Keplerian orbits by implementing a **Navarro-Frenk-White (NFW) Dark Matter Halo**. This allows for the simulation of realistic galaxies where invisible mass dictates the motion of outer stars, maintaining a flat rotation curve.

$$M(\lt r) = 4 \pi \rho_0 R_s^3 \left[ \ln\left(1 + \frac{r}{R_s}\right) - \frac{r/R_s}{1 + r/R_s} \right]$$

---

## Installation & Build

### Prerequisites
* **C++ Compiler** (GCC, Clang) supporting C++17
* **CMake** (Version 3.14+)
* **Python 3.11**
* **Libraries:** `pygame`, `PyOpenGL`, `numpy`, `pybind11`

### 1. Python Environment Setup
```bash
python3 -m venv venv
source venv/bin/activate
pip install pygame PyOpenGL PyOpenGL_accelerate numpy pybind11
```

### 2. Compile the Engine

For Linux/Standard Environments (via CMake):
```bash
mkdir build && cd build
cmake -DOpenMP_ROOT=$(brew --prefix libomp) ..
cmake --build .
```

For macOS / Apple Silicon (via Manual Script):
Due to Apple Clang and Homebrew ABI constraints, a dedicated build script is provided to ensure dynamic lookup and OpenMP linking.

```bash
chmod +x compile_manual.sh
./compile_manual.sh
```

### 3. Running the Simulation**

The engine is controlled via JSON configuration files, allowing rapid iteration of galactic parameters without recompiling.

```bash
python3 realtime_viz.py default_config.json
```

Interactive Controls

-- Left Click & Drag: Pan the camera

-- Scroll Wheel: Zoom in/out

-- 'T' Key: Toggle real-time Barnes-Hut Quadtree visualization grid

-- 'ESC': Exit simulation

