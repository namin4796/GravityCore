#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <cmath>
#include <memory>
#include <algorithm>
#include <iostream>

namespace py = pybind11;

// --- QuadTree Implementation ---
// The Tree Builder will be a helper class 
// that has access to the particle arrays.

class BarnesHutTree {
    struct Node {
        double x_min, y_min, size;
        double mass = 0;
        double cx = 0, cy = 0;
        bool is_leaf = true;

        int children[4] = {-1, -1, -1, -1};
        int particle_idx = -1; // -1 if empty or internal
        //std::unique_ptr<Node> children[4]; // NW, NE, SW, SE

        Node(double x, double y, double s) : x_min(x), y_min(y), size(s) {}
    };

    std::vector<Node> nodes;

    //std::unique_ptr<Node> root;
    const std::vector<double>& px;
    const std::vector<double>& py;
    const std::vector<double>& pm;

    std::vector<double>& debug_rects;

public:
    BarnesHutTree(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& m,
            std::vector<double>& rects_out) 
        : px(x), py(y), pm(m), debug_rects(rects_out) {

        debug_rects.clear(); // Reset for this frame
    
        // reserve memory to prevent reallocation
        nodes.reserve(px.size() * 4);    

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
        
        // create root (index 0)
        //root = std::make_unique<Node>(mid_x - size/2.0, mid_y - size/2.0, size);
        nodes.emplace_back(mid_x - size/2.0, mid_y - size/2.0, size);
            
        // 2. Insert All Particles
        for (size_t i = 0; i < px.size(); ++i) {
            insert(0, i); // begin at root (0)
        }
        // Collect boxes for visualization
        collect_boxes(0);
    }

    // Helper: create a child node and return its index
    int create_child(double x, double y, double s) {
        nodes.emplace_back(x, y, s);
        return nodes.size() - 1;
    }

    void collect_boxes(int node_idx) {
        if (node_idx == -1) return;

        const Node& node = nodes[node_idx];

        // Store format: [x, y, size]
        // // Only store if it has mass (is active)
        if (node.mass >0) {
            debug_rects.push_back(node.x_min);
            debug_rects.push_back(node.y_min);
            debug_rects.push_back(node.size);

            if (!node.is_leaf) {
                for(int child_idx : node.children) collect_boxes(child_idx);
            }
        }
    }

    void insert(int node_idx, int idx) {
        // we use reference wrapper or direct access carefully.
        // Warning: vector resizing invalidates references, so indices are used.

        // Update Center of Mass
        double p_mass = pm[idx];
        double p_x = px[idx];
        double p_y = py[idx];

        // access node by index
        Node& node = nodes[node_idx];

        if (node.mass == 0) {
            node.cx = p_x;
            node.cy = p_y;
            node.mass = p_mass;
            node.particle_idx = idx;
            return; // It was empty, now it's a leaf with this particle
        }

        // Standard COM update
        double total_m = node.mass + p_mass;
        node.cx = (node.cx * node.mass + p_x * p_mass) / total_m;
        node.cy = (node.cy * node.mass + p_y * p_mass) / total_m;
        node.mass = total_m;

        // If it's a leaf and contains a particle, we must split
        if (node.is_leaf) {
             // If we are landing exactly on top of another particle, skip to avoid infinite depth
            if (std::abs(px[node.particle_idx] - p_x) < 1e-7 &&
                std::abs(py[node.particle_idx] - p_y) < 1e-7) {
                 return; 
            }
            
            node.is_leaf = false;
            // Push the *existing* particle down
            int old_idx = node.particle_idx;
            node.particle_idx = -1; // It is no longer a leaf holding a particle
            push_to_child(node_idx, old_idx); // push old particle
        }

        // If not a leaf, push the *new* particle down
        push_to_child(node_idx, idx); // push new particle
    }

    void push_to_child(int node_idx, int idx) {
        double p_x = px[idx];
        double p_y = py[idx];
        double mid_x = nodes[node_idx].x_min + nodes[node_idx].size / 2.0;
        double mid_y = nodes[node_idx].y_min + nodes[node_idx].size / 2.0;

        int quad = 0;
        if (p_x >= mid_x) quad += 1; // East
        if (p_y >= mid_y) quad += 2; // South

        if (nodes[node_idx].children[quad] == -1) {
            double sub_size = nodes[node_idx].size / 2.0;
            double sub_x = (quad % 2 == 0) ? nodes[node_idx].x_min : mid_x;
            double sub_y = (quad < 2) ? nodes[node_idx].y_min : mid_y;

            //create child in the pool
            int child_idx = create_child(sub_x, sub_y, sub_size);

            nodes[node_idx].children[quad] = child_idx;
        }
        insert(nodes[node_idx].children[quad], idx);
    }

    void compute_force(int idx, double& fx, double& fy, double G, double theta) {
        calculate_force_recursive(0, idx, fx, fy, G, theta);
    }

    void calculate_force_recursive(int node_idx, int target_idx, double& fx, double& fy, double G, double theta) {
        if (node_idx == -1) return;


        const Node& node = nodes[node_idx];
        if(node.mass == 0) return;

        double dx = node.cx - px[target_idx];
        double dy = node.cy - py[target_idx];
        double dist_sq = dx*dx + dy*dy;
        double dist = std::sqrt(dist_sq);

        // Self-interaction check
        if (dist < 1e-9) return;

        // Barnes-Hut Criterion: s / d < theta
        // s = width of region, d = distance to center of mass
        bool far_enough = (node.size / dist) < theta;

        if (node.is_leaf || far_enough) {
            // Treat as single body
            double softening = 0.01; // Avoid singularities
            double eff_dist_sq = dist_sq + softening; 
            double eff_dist = std::sqrt(eff_dist_sq);
            
            double f = (G * pm[target_idx] * node.mass) / (eff_dist_sq);
            
            fx += f * (dx / eff_dist);
            fy += f * (dy / eff_dist);
        } else {
            // Too close, open the node
            for (int i = 0; i < 4; i++) {
                if (node.children[i]) {
                    calculate_force_recursive(node.children[i], target_idx, fx, fy, G, theta);
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
        Universe(int n_particles, double rs, double rho, double g_val, double time_step) :
            px(n_particles),  
            py(n_particles),
            vx(n_particles),
            vy(n_particles),
            ax(n_particles),
            ay(n_particles),
            mass(n_particles),
            R_SCALE(rs),
            RHO_0(rho),
            G(g_val),
            dt(time_step) {} 

        void set_state(const std::vector<double>& x, const
                std::vector<double>& y, const std::vector<double>& vx_in,
                const std::vector<double>& vy_in, const std::vector<double>& m) {
            
            px = x; py = y; vx = vx_in; vy = vy_in; mass = m;
            
            // Initial force calculation
            compute_forces();
        }

        std::vector<double> get_x() { return px; }
        std::vector<double> get_y() { return py; }
        std::vector<double> get_vx() { return vx; }
        std::vector<double> get_vy() { return vy; }
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
        .def("get_vx", &Universe::get_vx)
        .def("get_vy", &Universe::get_vy)
        .def("get_tree_rects", &Universe::get_tree_rects) // EXPOSE IT HERE
        .def("step", &Universe::step);
}
