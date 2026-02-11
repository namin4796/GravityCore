#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "gravity.hpp"

namespace py = pybind11;

PYBIND11_MODULE(gravity_core, m) {
    m.doc() = "High-performance N-body gravity engine with Barnes-Hut";
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
