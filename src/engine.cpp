#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "gravity.hpp"

namespace py = pybind11;

PYBIND11_MODULE(gravity_core, m) {
    m.doc() = "High-performance N-body gravity engine with Barnes-Hut";
    py::class_<Universe>(m, "Universe")
        .def(py::init<int, double, double, double, double>())
        //.def("set_state", &Universe::set_state)
        //.def("get_x", &Universe::get_x)
        //.def("get_y", &Universe::get_y)
        //.def("get_vx", &Universe::get_vx)
        //.def("get_vy", &Universe::get_vy)
        //.def("get_tree_rects", &Universe::get_tree_rects) // EXPOSE IT HERE
        .def("init_galaxy", &Universe::init_galaxy)
        .def("step", &Universe::step)
        .def("get_positions", [](Universe &self) {
                // to convert SoA to an Nx2 array for OpenGL
                auto x = self.get_x();
                auto y = self.get_y();
                auto result = py::array_t<double>({ (ssize_t)x.size(), (ssize_t)2 });
                auto buf = result.mutable_unchecked<2>();
                for (size_t i = 0; i < x.size(); i++) {
                    buf(i, 0) = x[i];
                    buf(i, 1) = y[i];
                }
                return result;
        });
}
