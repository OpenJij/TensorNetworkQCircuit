#include <circuit_topology.hpp>
#include <pybind11/pybind11.h>

namespace qcircuit {
    namespace py = pybind11;

    void init_circuit_topology(py::module& m) {
        py::class_<CircuitTopology>(m, "CircuitTopology")
            .def(py::init<size_t>())
            .def("generate_link", &CircuitTopology::generateLink)
            .def("number_of_bits", &CircuitTopology::numberOfBits)
            .def("number_of_links", &CircuitTopology::numberOfLinks)
            .def("convert_to_dot_string", &CircuitTopology::convertToDotString,
                 py::arg("layout") = "neato",
                 py::arg("shape") = "circle");
    }
}
