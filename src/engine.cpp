#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <cmath>
#include <memory>
#include <algorithm>
#include <iostream>

namespace py = pybind11;

// --- QuadTree Implementation ---
/*
struct QuadNode {
    double x_min, y_min, size;  // Spatial bounds
    double mass;
    double center_x, center_y;  // Center of Mass
    
    bool is_leaf;
    std::vector<int> particle_indices; // Indices of particles in this node (if leaf)
    
    // Children: NW, NE, SW, SE
    std::unique_ptr<QuadNode> children[4];

    QuadNode(double x, double y, double s) 
        : x_min(x), y_min(y), size(s), mass(0), center_x(0), center_y(0), is_leaf(true) {
            for(int i=0; i<4; i++) children[i] = nullptr;
    }

    // Determine which quadrant a point belongs to
    int get_quadrant(double x, double y) {
        double mid_x = x_min + size / 2.0;
        double mid_y = y_min + size / 2.0;
        
        int quad = 0;
        if (x >= mid_x) quad += 1; // East
        if (y >= mid_y) quad += 2; // South
        return quad;
    }

    void insert(int p_idx, double p_x, double p_y, double p_m) {
        // 1. Update Center of Mass for this node
        if (mass == 0) {
            center_x = p_x;
            center_y = p_y;
            mass = p_m;
        } else {
            double total_m = mass + p_m;
            center_x = (center_x * mass + p_x * p_m) / total_m;
            center_y = (center_y * mass + p_y * p_m) / total_m;
            mass = total_m;
        }

        // 2. Leaf Logic
        if (is_leaf) {
            // If this is the first particle, just store it
            if (particle_indices.empty()) {
                particle_indices.push_back(p_idx);
                return;
            }

            // If we already have a particle, we must subdivide (unless very close/same spot)
            // Prevent infinite recursion on overlapping particles
            double dx = std::abs(p_x - center_x);
            double dy = std::abs(p_y - center_y);
            if(dx < 1e-9 && dy < 1e-9) {
                 particle_indices.push_back(p_idx);
                 return;
            }

            is_leaf = false;
            
            // Push existing particles down to children
            for (int existing_idx : particle_indices) {
                // We need to retrieve the position of the existing particle? 
                // Optimization: We don't store pos in node, so we can't push down easily without 
                // passing the full position arrays.
                // FIX: For this simplified version, let's just assume we only split on collision.
                // To do this strictly, we need access to the global arrays here.
                // See "Note on Implementation" below.
            }
            // For a robust implementation without global access in struct, 
            // we usually limit leaf capacity to 1 OR pass arrays. 
            // Let's use a simpler recursive approach in the Universe class helper.
        }
    }
;
*/
// We will use a slightly different approach: The Tree Builder will be a helper class 
// that has access to the particle arrays.

class BarnesHutTree {
    struct Node {
        double x_min, y_min, size;
        double mass = 0;
        double cx = 0, cy = 0;
        bool is_leaf = true;
        int particle_idx = -1; // -1 if empty or internal
        std::unique_ptr<Node> children[4]; // NW, NE, SW, SE

        Node(double x, double y, double s) : x_min(x), y_min(y), size(s) {}
    };

    std::unique_ptr<Node> root;
    const std::vector<double>& px;
    const std::vector<double>& py;
    const std::vector<double>& pm;

    std::vector<double>& debug_rects;

public:
    BarnesHutTree(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& m,
            std::vector<double>& rects_out) 
        : px(x), py(y), pm(m), debug_rects(rects_out) {

        debug_rects.clear(); // Reset for this frame
        
        // 1. Find Bounds
        double min_x = *std::min_element(px.begin(), px.end());
        double max_x = *std::max_element(px.begin(), px.end());
        double min_y = *std::min_element(py.begin(), py.end());
        double max_y = *std::max_element(py.begin(), py.end());

        double size_x = max_x - min_x;
        double size_y = max_y - min_y;
        double size = std::max(size_x, size_y) * 1.01; // slightly larger to avoid boundary edge cases
        
        // Center the square
        double mid_x = (min_x + max_x) / 2.0;
        double mid_y = (min_y + max_y) / 2.0;
        
        root = std::make_unique<Node>(mid_x - size/2.0, mid_y - size/2.0, size);

        // 2. Insert All Particles
        for (size_t i = 0; i < px.size(); ++i) {
            insert(root.get(), i);
        }
        // Collect boxes for visualization
        collect_boxes(root.get());
    }

    void collect_boxes(Node* node) {
        if (!node) return;

        // Store format: [x, y, size]
        // // Only store if it has mass (is active)
        if (node->mass >0) {
            debug_rects.push_back(node->x_min);
            debug_rects.push_back(node->y_min);
            debug_rects.push_back(node->size);

            if (!node->is_leaf) {
                for(auto& child : node->children) collect_boxes(child.get());
            }
        }
    }

    void insert(Node* node, int idx) {
        // Update Center of Mass
        double p_mass = pm[idx];
        double p_x = px[idx];
        double p_y = py[idx];

        if (node->mass == 0) {
            node->cx = p_x;
            node->cy = p_y;
            node->mass = p_mass;
            node->particle_idx = idx;
            return; // It was empty, now it's a leaf with this particle
        }

        // Standard COM update
        double total_m = node->mass + p_mass;
        node->cx = (node->cx * node->mass + p_x * p_mass) / total_m;
        node->cy = (node->cy * node->mass + p_y * p_mass) / total_m;
        node->mass = total_m;

        // If it's a leaf and contains a particle, we must split
        if (node->is_leaf) {
             // If we are landing exactly on top of another particle, skip to avoid infinite depth
            if (std::abs(px[node->particle_idx] - p_x) < 1e-7 && std::abs(py[node->particle_idx] - p_y) < 1e-7) {
                 return; 
            }
            
            node->is_leaf = false;
            // Push the *existing* particle down
            int old_idx = node->particle_idx;
            node->particle_idx = -1; // It is no longer a leaf holding a particle
            push_to_child(node, old_idx);
        }

        // If not a leaf, push the *new* particle down
        push_to_child(node, idx);
    }

    void push_to_child(Node* node, int idx) {
        double p_x = px[idx];
        double p_y = py[idx];
        double mid_x = node->x_min + node->size / 2.0;
        double mid_y = node->y_min + node->size / 2.0;

        int quad = 0;
        if (p_x >= mid_x) quad += 1; // East
        if (p_y >= mid_y) quad += 2; // South

        if (!node->children[quad]) {
            double sub_size = node->size / 2.0;
            double sub_x = (quad % 2 == 0) ? node->x_min : mid_x;
            double sub_y = (quad < 2) ? node->y_min : mid_y;
            node->children[quad] = std::make_unique<Node>(sub_x, sub_y, sub_size);
        }
        insert(node->children[quad].get(), idx);
    }

    void compute_force(int idx, double& fx, double& fy, double G, double theta) {
        calculate_force_recursive(root.get(), idx, fx, fy, G, theta);
    }

    void calculate_force_recursive(Node* node, int target_idx, double& fx, double& fy, double G, double theta) {
        if (!node || node->mass == 0) return;

        double dx = node->cx - px[target_idx];
        double dy = node->cy - py[target_idx];
        double dist_sq = dx*dx + dy*dy;
        double dist = std::sqrt(dist_sq);

        // Self-interaction check
        if (dist < 1e-9) return;

        // Barnes-Hut Criterion: s / d < theta
        // s = width of region, d = distance to center of mass
        bool far_enough = (node->size / dist) < theta;

        if (node->is_leaf || far_enough) {
            // Treat as single body
            double softening = 0.01; // Avoid singularities
            double eff_dist_sq = dist_sq + softening; 
            double eff_dist = std::sqrt(eff_dist_sq);
            
            double f = (G * pm[target_idx] * node->mass) / (eff_dist_sq);
            
            fx += f * (dx / eff_dist);
            fy += f * (dy / eff_dist);
        } else {
            // Too close, open the node
            for (int i = 0; i < 4; i++) {
                if (node->children[i]) {
                    calculate_force_recursive(node->children[i].get(), target_idx, fx, fy, G, theta);
                }
            }
        }
    }
};

// --- Main Engine ---

class Universe {
    std::vector<double> px, py; 
    std::vector<double> vx, vy; 
    std::vector<double> ax, ay; 
    std::vector<double> mass;
    std::vector<double> current_rects;
    double dt;
    double R_SCALE, RHO_0, G;

    // Tuning parameter for Barnes-Hut
    double THETA = 0.5; 

    private:
       void compute_forces() {
           size_t n = px.size();

           std::fill(ax.begin(), ax.end(), 0.0);
           std::fill(ay.begin(), ay.end(), 0.0);

           // 1. Dark Matter Halo Force (Analytic - stays the same)
           double M_factor = 4.0 * M_PI * RHO_0 * std::pow(R_SCALE, 3);

           #pragma omp parallel for
           for (size_t i = 0; i < n; i++) {
               double dx = px[i];
               double dy = py[i];
               double r_sq = dx*dx + dy*dy;
               double r = std::sqrt(r_sq);

               if (r > 1e-5) {
                   double x = r / R_SCALE;
                   double m_NFW = M_factor * (std::log(1.0+x) - (x/(1.0+x)));
                   double f_NFW = G * (m_NFW / r_sq);
                   ax[i] -= f_NFW * (dx/r);
                   ay[i] -= f_NFW * (dy/r);
               }
           }

           // 2. N-Body Gravity (NOW USING BARNES-HUT)
           // Build the tree
           BarnesHutTree tree(px, py, mass, current_rects);

           // Calculate forces using the tree
           // Note: Tree building is serial, but force calculation can be parallelized easily!
           #pragma omp parallel for
           for (size_t i = 0; i < n; i++) {
               double fx = 0.0;
               double fy = 0.0;
               
               tree.compute_force(i, fx, fy, G, THETA);

               ax[i] += fx / mass[i];
               ay[i] += fy / mass[i];
           }
       }

    public:
        Universe(int n_particles, double rs, double rho, double g_val, double time_step) {
            R_SCALE = rs;
            RHO_0 = rho;
            G = g_val;
            dt = time_step;

            px.resize(n_particles, 0.0);
            py.resize(n_particles, 0.0);
            vx.resize(n_particles, 0.0);
            vy.resize(n_particles, 0.0);
            ax.resize(n_particles, 0.0);
            ay.resize(n_particles, 0.0);
            mass.resize(n_particles, 1.0);
        }

        void set_state(const std::vector<double>& x, const
                std::vector<double>& y, const std::vector<double>& vx_in,
                const std::vector<double>& vy_in, const std::vector<double>& m) {
            
            if (x.size() != m.size()) {
                throw std::runtime_error("Input array sizes do not match!");
            }
            px = x; py = y; vx = vx_in; vy = vy_in; mass = m;
            ax.resize(px.size(), 0.0);
            ay.resize(py.size(), 0.0);
            
            // Initial force calculation
            compute_forces();
        }

        std::vector<double> get_x() { return px; }
        std::vector<double> get_y() { return py; }
//        std::vector<double> get_vx() { return vx; }
//        std::vector<double> get_vy() { return vy; }
        std::vector<double> get_tree_rects() { return current_rects; }

       void step() {
           size_t n = px.size();
           // Velocity Verlet: First half-kick and drift
           #pragma omp parallel for
           for (size_t i = 1; i < n; i++) {
                px[i] += vx[i] * dt + 0.5 * ax[i] * dt * dt;         
                py[i] += vy[i] * dt + 0.5 * ay[i] * dt * dt;         
                vx[i] += 0.5 * ax[i] * dt;
                vy[i] += 0.5 * ay[i] * dt;
           }

           // Recompute forces with new positions
           compute_forces();

           // Final half-kick
           #pragma omp parallel for
           for (size_t i = 1; i < n; i++) {
               vx[i] += 0.5 * ax[i] * dt;
               vy[i] += 0.5 * ay[i] * dt;
           }
       }
};

PYBIND11_MODULE(gravity_core, m) {
    m.doc() = "high-performance N-body gravity engine with Barnes-Hut";
    py::class_<Universe>(m, "Universe")
        .def(py::init<int, double, double, double, double>())
        .def("set_state", &Universe::set_state)
        .def("get_x", &Universe::get_x)
        .def("get_y", &Universe::get_y)
        //.def("get_vx", &Universe::get_vx)
        //.def("get_vy", &Universe::get_vy)
        .def("get_tree_rects", &Universe::get_tree_rects) // EXPOSE IT HERE
        .def("step", &Universe::step);
}
