# GravityCore: High-Performance N-Body Physics Engine

![Language](https://img.shields.io/badge/language-C%2B%2B17%20%7C%20Python-blue)
![Build](https://img.shields.io/badge/build-CMake-green)
![License](https://img.shields.io/badge/license-MIT-lightgrey)

**GravityCore** is a hybrid physics simulation engine capable of calculating gravitational interactions for thousands of bodies in real-time. It leverages **Modern C++ (C++17)** for high-performance computing and **Python** for high-level orchestration and visualization.

The project demonstrates the power of **Data-Oriented Design**, achieving a **~30x speedup** over pure Python implementations by optimizing memory access patterns and utilizing hardware parallelism.

![Rotation Curve](https://github.com/namin4796/GravityCore/blob/main/rotation_curve.png)
* Figure: Comparison of rotation curve of stars with and without the NFW dark matter.*
---

## 🚀 Key Features

### 1. High-Performance Architecture
* **Hybrid C++/Python:** Computationally expensive $O(N^2)$ force calculations are offloaded to a compiled C++ backend, while simulation setup and analysis remain in Python.
* **Parallel Computing:** Utilizes **OpenMP** to parallelize force loops across all available CPU cores.
* **Data-Oriented Design (SoA):** Implements a **Structure of Arrays** memory layout (instead of Array of Structures) to maximize CPU cache locality and vectorization potential.

### 2. Scientific Accuracy
* **Symplectic Integration:** Uses the **Velocity Verlet** integrator to ensure energy conservation over long simulation times (unlike standard Euler integration).
* **Robust Physics:** Includes gravitational softening parameters ($\epsilon$) to prevent numerical singularities during close-range stellar collisions.
* **N-Body Interaction:** Simulates full all-to-all gravitational forces, allowing for complex phenomena like star clustering, tidal stripping, and galaxy formation.

---

## 🛠️ Installation & Build

### Prerequisites
* **C++ Compiler** (GCC, Clang, or MSVC) supporting C++17
* **CMake** (Version 3.14+)
* **Python 3.x**
* **Libraries:** `numpy`, `matplotlib` (for visualization)

### Build Instructions
The project uses **CMake** to automatically fetch dependencies (like `pybind11`) and compile the Python module.

```bash
# 1. Clone the repository
git clone [https://github.com/YOUR_USERNAME/GravityCore.git](https://github.com/YOUR_USERNAME/GravityCore.git)
cd GravityCore

# 2. Create build directory
mkdir build && cd build

# 3. Configure and Compile
cmake ..
cmake --build .
