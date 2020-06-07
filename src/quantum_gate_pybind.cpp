#include <quantum_gate.hpp>
#include <pybind11/pybind11.h>

namespace qcircuit {
    namespace py = pybind11;

    void init_quantum_gate(py::module& m) {
        py::class_<OneSiteGate>(m, "OneSiteGate");
        py::class_<TwoSiteGate>(m, "TwoSiteGate");

        py::class_<Id, OneSiteGate>(m, "Id").def(py::init<size_t>());
        py::class_<X, OneSiteGate>(m, "X").def(py::init<size_t>());
        py::class_<Y, OneSiteGate>(m, "Y").def(py::init<size_t>());
        py::class_<Z, OneSiteGate>(m, "Z").def(py::init<size_t>());
        py::class_<Proj_0, OneSiteGate>(m, "Proj_0").def(py::init<size_t>());
        py::class_<Proj_1, OneSiteGate>(m, "Proj_1").def(py::init<size_t>());
        py::class_<Proj_0_to_1, OneSiteGate>(m, "Proj_0_to_1").def(py::init<size_t>());
        py::class_<Proj_1_to_0, OneSiteGate>(m, "Proj_1_to_0").def(py::init<size_t>());
        py::class_<H, OneSiteGate>(m, "H").def(py::init<size_t>());

        py::class_<CNOT, TwoSiteGate>(m, "CNOT").def(py::init<size_t, size_t>());
        py::class_<CY, TwoSiteGate>(m, "CY").def(py::init<size_t, size_t>());
        py::class_<CZ, TwoSiteGate>(m, "CZ").def(py::init<size_t, size_t>());
        py::class_<Swap, TwoSiteGate>(m, "Swap").def(py::init<size_t, size_t>());
    }
}
