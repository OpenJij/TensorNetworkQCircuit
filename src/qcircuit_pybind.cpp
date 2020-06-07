#include <itensor/all.h>
#include <qcircuit.hpp>
#include <circuit_topology.hpp>
#include <quantum_gate.hpp>
#include <pybind11/pybind11.h>

namespace qcircuit {
    namespace py = pybind11;

    void init_qcircuit(py::module& m) {
        py::class_<QCircuit>(m, "QCircuit")
            .def(py::init<const CircuitTopology&>())
            .def("apply", [](QCircuit& self, const OneSiteGate& g1, const OneSiteGate& g2) {
                              self.apply(g1, g2);
                          })
            .def("apply", [](QCircuit& self, const TwoSiteGate& g) {
                              self.apply(g);
                          })
            .def("probabilityOfZero", &QCircuit::probabilityOfZero)
            .def("observeQubit", [](QCircuit& self, size_t ind) {
                                     return self.observeQubit(ind);
                                 });
        /*
         * In above code, some member functions are wrapped in lambdas.
         * This is because pybind can't bind functions with default arguments correctly
         * (There might be better solution).
         */
    }
}
