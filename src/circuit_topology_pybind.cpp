#include <circuit_topology.hpp>
#include <pybind11/pybind11.h>

namespace qcircuit {
    namespace py = pybind11;

    void init_circuit_topology(py::module& m) {
        py::class_<CircuitTopology>(m, "CircuitTopology")
            .def(py::init<size_t>())
            .def("generateLink", &CircuitTopology::generateLink)
            .def("numberOfBits", &CircuitTopology::numberOfBits)
            .def("numberOfLinks", &CircuitTopology::numberOfLinks)
            .def("exportDotString", &CircuitTopology::exportDotString,
                 py::arg("filename"),
                 py::arg("layout") = "neato",
                 py::arg("shape") = "circle");
    }
}
