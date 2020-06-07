#include <pybind11/pybind11.h>

namespace qcircuit {
    namespace py = pybind11;

    void init_qcircuit(py::module&);
    void init_circuit_topology(py::module&);
    void init_circuits(py::module&);
    void init_quantum_gate(py::module&);

    PYBIND11_MODULE(qcircuit, m) {
        init_qcircuit(m);
        init_circuit_topology(m);
        init_circuits(m);
        init_quantum_gate(m);
    }
}
