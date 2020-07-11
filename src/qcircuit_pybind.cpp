#include <vector>
#include <itensor/all.h>
#include <qcircuit.hpp>
#include <circuit_topology.hpp>
#include <quantum_gate.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace qcircuit {
    namespace py = pybind11;

    void init_qcircuit(py::module& m) {
        py::class_<QCircuit>(m, "QCircuit")
            .def(py::init<const CircuitTopology&>())
            .def("apply", py::overload_cast<const OneSiteGate&, const OneSiteGate&>(&QCircuit::apply))
            .def("apply", py::overload_cast<const OneSiteGate&>(&QCircuit::apply))
            .def("apply", py::overload_cast<const TwoSiteGate&>(&QCircuit::apply))
            .def("get_cursor", &QCircuit::getCursor)
            .def("move_cursor_along", py::overload_cast<const std::vector<size_t>&>(&QCircuit::moveCursorAlong))
            .def("probability_of_zero", &QCircuit::probabilityOfZero)
            .def("observe_qubit", py::overload_cast<size_t>(&QCircuit::observeQubit))
            .def_property("cutoff", &QCircuit::getCutoff, &QCircuit::setCutoff)
            .def_property("max_dim", &QCircuit::getMaxDim, &QCircuit::setMaxDim);
    }
}
