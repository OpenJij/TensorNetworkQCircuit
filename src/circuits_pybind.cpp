#include <circuits.hpp>
#include <pybind11/pybind11.h>

namespace qcircuit {
    namespace py = pybind11;

    void init_circuits(py::module& m) {
        m.def("make_ibmq_topology", &make_ibmq_topology);
        m.def("make_chain", &make_chain,
              py::arg("size"),
              py::arg("periodic") = true);

        m.def("make_alltoall_topology", &make_alltoall_topology);
    }
}
