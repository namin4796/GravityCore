#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <cmath>

namespace py = pybind11;

class Universe {
    std::vector<double> px, py; // positions
    std::vector<double> vx, vy; // velocities
    std::vector<double> mass; // masses
    const double G = 1.0;//6.67430e-11;
    double dt = 0.005;

    public:
        Universe(int n_particles) {
            px.resize(n_particles, 0.0);
            px.resize(n_particles, 0.0);
            vx.resize(n_particles, 0.0);
            vx.resize(n_particles, 0.0);
            mass.resize(n_particles, 1.0);
        }

        // Init with data from Python
        void set_state(const std::vector<double>& x, const
                std::vector<double>& y, const std::vector<double>& vx_in,
                const std::vector<double>& vy_in, const std::vector<double>& m) {

            if (x.size() != y.size() || x.size() != vx_in.size() || x.size() != vy_in.size()
                    || x.size() != m.size()) {
                throw std::runtime_error("Input array sizes (x, y, vx_in, vy_in, mass) do not match!");
            }
            px = x;
            py = y;
            vx = vx_in;
            vy = vy_in;
            mass = m;

            //resize velocities
            vx.resize(px.size(), 0.0);
            vy.resize(px.size(), 0.0);

            //vx.assign(px.size(), 0.1);
        }

        // Get back data from Python for plotting
        std::vector<double> get_x() { return px; }
        std::vector<double> get_y() { return py; }
                           
       // Grav. force calculation O(N^2) 
       void step() {
           size_t n = px.size();

           //calculate forces / update velocities
           #pragma omp parallel for
           for (size_t i = 0; i < n; i++) {
               double fx = 0.0;
               double fy = 0.0;

               for (size_t j = 0; j < n; j++) {
                   if (i == j)
                       continue;

                   double dx = px[j] - px[i];
                   double dy = py[j] - py[i];
                   double dist_sq = dx*dx + dy*dy + 1e-9; // +epsilon to avoid div by 0
                   double dist = std::sqrt(dist_sq);

                   // grav for F = G*(m1*m2)/r2
                   double f = (G * mass[i] * mass[j]) / dist_sq;

                   fx += f * (dx / dist);
                   fy += f * (dy / dist);
               }

               //F = ma -> a = F/m
               vx[i] += (fx / mass[i]) * dt;
               vy[i] += (fy / mass[i]) * dt;
           }

           // update positions
           for (size_t i = 0; i < n; i++) {
               px[i] += vx[i] * dt;
               py[i] += vy[i] * dt;
           }
       }
};

//Bind it to Python
PYBIND11_MODULE(gravity_core, m) {
    m.doc() = "high-performance N-body gravity engine";

    py::class_<Universe>(m, "Universe")
        .def(py::init<int>())
        .def("set_state", &Universe::set_state)
        .def("get_x", &Universe::get_x)
        .def("get_y", &Universe::get_y)
        .def("step", &Universe::step);
}

