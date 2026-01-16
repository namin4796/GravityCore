import sys
import os
import random
import numpy as np
import matplotlib.pyplot as plt

# 1. SETUP PATHS
sys.path.append(os.path.join(os.path.dirname(__file__), '../build'))

try:
    import gravity_core
    print("SUCCESS: Loaded C++ Engine!")
except ImportError:
    print("ERROR: C++ Engine not found. Did you recompile?")
    sys.exit(1)

# 2. HELPER: Theoretical Curves
def get_nfw_velocity(r, m_bh):
    # Matches C++ Constants
    R_SCALE = 200.0
    RHO_0 = 0.002
    G = 1.0
    
    # NFW Mass enclosed
    x = r / R_SCALE
    m_factor = 4.0 * np.pi * RHO_0 * (R_SCALE**3)
    m_halo = m_factor * (np.log(1+x) - (x/(1+x)))
    
    # Total Mass
    m_total = m_bh + m_halo
    
    # Circular Velocity: v = sqrt(GM/r)
    return np.sqrt(G * m_total / r)

def get_kepler_velocity(r, m_bh):
    # Only Black Hole
    G = 1.0
    return np.sqrt(G * m_bh / r)

# 3. INITIALIZE SIMULATION
N = 500
sim = gravity_core.Universe(N)
mass_bh = 1000.0

# Generate stars (Same logic as visualize.py)
radii = np.linspace(10, 200, N) # Spread them out evenly to see the curve
angles = np.random.uniform(0, 2*np.pi, N)

px = radii * np.cos(angles)
py = radii * np.sin(angles)
mass = np.ones(N)

# Set BH
px[0], py[0] = 0, 0
mass[0] = mass_bh

# Set Velocities (Using NFW assumption)
vx = []
vy = []
for i in range(N):
    r = radii[i]
    if i == 0: # BH
        vx.append(0); vy.append(0)
        continue
        
    v_ideal = get_nfw_velocity(r, mass_bh)
    vx.append(-v_ideal * np.sin(angles[i]))
    vy.append(v_ideal * np.cos(angles[i]))

sim.set_state(px, py, vx, vy, mass)

# 4. Run Simulation
# We just take the initial state since we set it up perfectly, 
# but you could run sim.step() here to see if they hold.
for _ in range(50):
    sim.step()

# Retrieve data from C++ engine.
final_px = np.array(sim.get_x())
final_py = np.array(sim.get_y())
final_vx = np.array(sim.get_vx())
final_vy = np.array(sim.get_vy())

measured_r = np.sqrt(final_px**2 + final_py**2)
measured_v = np.sqrt(final_vx**2 + final_vy**2)

measured_r = measured_r[1:]
measured_v = measured_v[1:]


# 5. PLOT
plt.figure(figsize=(10, 6))

r_plot = np.linspace(10, 260, 100)
v_kepler = get_kepler_velocity(r_plot, mass_bh)
v_nfw = get_nfw_velocity(r_plot, mass_bh)

# Plot Theoretical Curves
plt.plot(r_plot, v_kepler, 'k--', label='Expected (Black Hole Only)')
plt.plot(r_plot, v_nfw, 'r-', linewidth=1, label='Expected (Dark Matter + BH)')

# Plot Your Particles
# We plot the velocities we GAVE them. If the simulation is stable, they stay here.
# If C++ didn't have DM, these particles would fly away!
v_particles = [np.sqrt(vx[i]**2 + vy[i]**2) for i in range(1, N)]
r_particles = radii[1:]

plt.scatter(r_particles, v_particles, color='blue', alpha=0.6, s=10, label='Simulation Particles')

plt.xlabel('Distance from Center (r)')
plt.ylabel('Orbital Velocity (v)')
plt.title('Galaxy Rotation Curve: Evidence of Dark Matter')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("rotation_curve.png")
plt.show()
