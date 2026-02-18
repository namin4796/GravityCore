#pragma once
#include <vector>
#include <cmath>
#include <random>
#include <memory>
#include <algorithm>
#include <iostream>
#include <omp.h>

// --- QuadTree Implementation (Linear Memory Pool) ---
class BarnesHutTree {
    struct Node {
        double x_min, y_min, size;
        double mass = 0;
        double cx = 0, cy = 0;
        bool is_leaf = true;
        int children[4] = {-1, -1, -1, -1}; 
        int particle_idx = -1; 
        Node(double x, double y, double s) : x_min(x), y_min(y), size(s) {}
    };

    std::vector<Node> nodes;
    const std::vector<double>& px;
    const std::vector<double>& py;
    const std::vector<double>& pm;
    std::vector<double>& debug_rects;

public:
    BarnesHutTree(const std::vector<double>& x, const std::vector<double>& y, 
                  const std::vector<double>& m, std::vector<double>& rects_out) 
        : px(x), py(y), pm(m), debug_rects(rects_out) {
        
        debug_rects.clear(); 
        nodes.reserve(px.size() * 4); 

        double min_x = *std::min_element(px.begin(), px.end());
        double max_x = *std::max_element(px.begin(), px.end());
        double min_y = *std::min_element(py.begin(), py.end());
        double max_y = *std::max_element(py.begin(), py.end());

        double size = std::max(max_x - min_x, max_y - min_y) * 1.01;
        double mid_x = (min_x + max_x) / 2.0;
        double mid_y = (min_y + max_y) / 2.0;

        nodes.emplace_back(mid_x - size/2.0, mid_y - size/2.0, size);

        for (size_t i = 0; i < px.size(); ++i) insert(0, i);
        collect_boxes(0);
    }

    int create_child(double x, double y, double s) {
        nodes.emplace_back(x, y, s);
        return nodes.size() - 1;
    }

    // FIX: Using indices instead of references to prevent invalidation during resize
    void insert(int node_idx, int p_idx) {
        double p_mass = pm[p_idx];
        double p_x = px[p_idx];
        double p_y = py[p_idx];

        // 1. Update Mass (Access via index)
        if (nodes[node_idx].mass == 0) {
            nodes[node_idx].cx = p_x; 
            nodes[node_idx].cy = p_y; 
            nodes[node_idx].mass = p_mass;
            nodes[node_idx].particle_idx = p_idx;
            return;
        }

        double total_m = nodes[node_idx].mass + p_mass;
        nodes[node_idx].cx = (nodes[node_idx].cx * nodes[node_idx].mass + p_x * p_mass) / total_m;
        nodes[node_idx].cy = (nodes[node_idx].cy * nodes[node_idx].mass + p_y * p_mass) / total_m;
        nodes[node_idx].mass = total_m;

        // 2. Split if Leaf
        if (nodes[node_idx].is_leaf) {
            // Prevent infinite recursion on overlapping particles
            int old_p = nodes[node_idx].particle_idx;
            if (std::abs(px[old_p] - p_x) < 1e-7 && std::abs(py[old_p] - p_y) < 1e-7) return;

            nodes[node_idx].is_leaf = false;
            nodes[node_idx].particle_idx = -1;
            
            // Push the old particle down (Recursive call might resize vector!)
            push_to_child(node_idx, old_p);
        }

        // Push the new particle down
        push_to_child(node_idx, p_idx); 
    }

    void push_to_child(int node_idx, int p_idx) {
        // Always re-read bounds from the vector, do not cache 'Node&'
        double mid_x = nodes[node_idx].x_min + nodes[node_idx].size / 2.0;
        double mid_y = nodes[node_idx].y_min + nodes[node_idx].size / 2.0;
        double p_x = px[p_idx];
        double p_y = py[p_idx];

        int quad = 0;
        if (p_x >= mid_x) quad += 1;
        if (p_y >= mid_y) quad += 2;

        if (nodes[node_idx].children[quad] == -1) {
            double sub_size = nodes[node_idx].size / 2.0;
            double sub_x = (quad % 2 == 0) ? nodes[node_idx].x_min : mid_x;
            double sub_y = (quad < 2) ? nodes[node_idx].y_min : mid_y;
            
            int child_idx = create_child(sub_x, sub_y, sub_size);
            
            // Re-access node_idx because create_child might have moved memory
            nodes[node_idx].children[quad] = child_idx;
        }

        insert(nodes[node_idx].children[quad], p_idx);
    }

    void collect_boxes(int node_idx) {
        if (node_idx == -1) return;
        // const ref is okay here because we are not adding/resizing
        const Node& node = nodes[node_idx];
        
        if (node.mass > 0) {
            debug_rects.push_back(node.x_min);
            debug_rects.push_back(node.y_min);
            debug_rects.push_back(node.size);
            
            if (!node.is_leaf) {
                for(int child_idx : node.children) collect_boxes(child_idx);
            }
        }
    }

    void compute_force(int p_idx, double& fx, double& fy, double G, double theta) {
        calculate_force_recursive(0, p_idx, fx, fy, G, theta);
    }

    void calculate_force_recursive(int node_idx, int target_idx, double& fx, double& fy, double G, double theta) {
        if (node_idx == -1) return;
        const Node& node = nodes[node_idx]; // Safe (read-only)
        
        if (node.mass == 0) return;

        double dx = node.cx - px[target_idx];
        double dy = node.cy - py[target_idx];
        double dist_sq = dx*dx + dy*dy;
        double dist = std::sqrt(dist_sq);

        if (dist < 1e-9) return;

        bool far_enough = (node.size / dist) < theta;

        if (node.is_leaf || far_enough) {
            double softening = 0.01;
            double eff_dist_sq = dist_sq + softening; 
            double eff_dist = std::sqrt(eff_dist_sq);
            double f = (G * pm[target_idx] * node.mass) / (eff_dist_sq);
            
            fx += f * (dx / eff_dist);
            fy += f * (dy / eff_dist);
        } else {
            for (int i = 0; i < 4; i++) {
                calculate_force_recursive(node.children[i], target_idx, fx, fy, G, theta);
            }
        }
    }
};

// --- Universe Class (Unchanged) ---
class Universe {
    std::vector<double> px, py, vx, vy, ax, ay, mass, current_rects;
    double dt, R_SCALE, RHO_0, G;
    double THETA = 0.5;

    void compute_forces() {
        size_t n = px.size();
        std::fill(ax.begin(), ax.end(), 0.0);
        std::fill(ay.begin(), ay.end(), 0.0);

        double M_factor = 4.0 * M_PI * RHO_0 * std::pow(R_SCALE, 3);
        #pragma omp parallel for
        for (size_t i = 0; i < n; i++) {
            double dx = px[i], dy = py[i];
            double r = std::sqrt(dx*dx + dy*dy);
            if (r > 1e-5) {
                double x = r / R_SCALE;
                double m_NFW = M_factor * (std::log(1.0+x) - (x/(1.0+x)));
                double f_NFW = G * (m_NFW / (r*r));
                ax[i] -= f_NFW * (dx/r);
                ay[i] -= f_NFW * (dy/r);
            }
        }

        BarnesHutTree tree(px, py, mass, current_rects);

        #pragma omp parallel for
        for (size_t i = 0; i < n; i++) {
            double fx = 0.0, fy = 0.0;
            tree.compute_force(i, fx, fy, G, THETA);
            ax[i] += fx / mass[i];
            ay[i] += fy / mass[i];
        }
    }

public:
    Universe(int n, double rs, double rho, double g, double t) 
        : px(n), py(n), vx(n), vy(n), ax(n), ay(n), mass(n), dt(t), R_SCALE(rs), RHO_0(rho), G(g) {}

    void set_state(const std::vector<double>& x, const std::vector<double>& y, 
                   const std::vector<double>& vx_in, const std::vector<double>& vy_in, 
                   const std::vector<double>& m) {
        px = x; py = y; vx = vx_in; vy = vy_in; mass = m;
        compute_forces();
    }
    
    std::vector<double> get_x() { return px; }
    std::vector<double> get_y() { return py; }
    std::vector<double> get_vx() { return vx; }
    std::vector<double> get_vy() { return vy; }
    std::vector<double> get_mass() { return mass; }
    std::vector<double> get_tree_rects() { return current_rects; }

    void step() {
        size_t n = px.size();
        #pragma omp parallel for
        for (size_t i = 1; i < n; i++) {
            px[i] += vx[i]*dt + 0.5*ax[i]*dt*dt;
            py[i] += vy[i]*dt + 0.5*ay[i]*dt*dt;
            vx[i] += 0.5*ax[i]*dt;
            vy[i] += 0.5*ay[i]*dt;
        }
        compute_forces();
        #pragma omp parallel for
        for (size_t i = 1; i < n; i++) {
            vx[i] += 0.5*ax[i]*dt;
            vy[i] += 0.5*ay[i]*dt;
        }

        //to account for absorption of the stars by BH
        handle_collisions();
    }

    void init_galaxy(double radius, double central_mass) {
        std::mt19937 gen(42);
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        px[0] = 0;
        py[0] = 0;
        vx[0] = 0;
        vy[0] = 0;
        mass[0] = central_mass;

        for (size_t i = 1; i < px.size(); ++i) {
            double r = radius * (0.1 + 0.9*std::sqrt(dist(gen)));
            double theta = dist(gen) * 2.0 * M_PI;

            px[i] = r*std::cos(theta);
            py[i] = r*std::sin(theta);
            mass[i] = 1.0 + dist(gen)*2.0;

            double v_mag = std::sqrt(G * central_mass / r);
            vx[i] = -v_mag * std::sin(theta);
            vy[i] = v_mag * std::cos(theta);
        }
    }

    void handle_collisions() {
        double eh_radius = 5.0; // Event Horizon of the central BH
        
        for (size_t i = 1; i < px.size(); ) {
            double dx = px[i] - px[0];
            double dy = py[i] - py[0];
            double dist_sq = dx*dx + dy*dy;

            if (dist_sq < (eh_radius * eh_radius)) {
                // if star comes close to the BH, its mass increases
                // and acquires velocity for conserving momentum :
                // v_new = (m1*v1 + m2*v2) / (m1 + m2)
                double new_mass = mass[0] + mass[i];
                vx[0] = (mass[0]*vx[0] + mass[i]*vx[i]) / new_mass;
                vy[0] = (mass[0]*vy[0] + mass[i]*vy[i]) / new_mass;
                mass[0] = new_mass;

                // delete the star that is absorbed by the BH
                size_t last = px.size() - 1;
                if (i != last) {
                    px[i] = px[last];
                    py[i] = py[last];
                    vx[i] = vx[last];
                    vy[i] = vy[last];
                    ax[i] = ax[last];
                    ay[i] = ay[last];
                    mass[i] = mass[last];
                }

                //remoe the last element to avoid duplication
                px.pop_back();
                py.pop_back();
                vx.pop_back();
                vy.pop_back();
                ax.pop_back();
                ay.pop_back();
                mass.pop_back();

            }
            else {
                i++;
            }
        }
    }                                     
};
