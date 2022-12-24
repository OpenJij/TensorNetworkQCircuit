# README for developers

## Dependencies
[ITensor](https://itensor.org/) v3 and [pybind11](https://github.com/pybind/pybind11)
are required.

Directory structure which is not contained in Git repository
is currently as follows:
```
.
├── build (Working directory for CMake out-of-source build)
├── html (Doxygen output)
├── external
    ├── itensor
    └── pybind11
```

## Install

### From GitHub
```sh
mkdir external
cd external
git clone https://github.com/ITensor/ITensor itensor
cd itensor
# configure `options.mk` and `make` following the ITensor install document https://github.com/ITensor/ITensor/blob/v3/INSTALL.md
cd ..
# Install on your virtual environment
pip install .
# or create wheel package
pip install wheel
python setup.py bdist_wheel
```

## C++ Example

```c++
#include <itensor/all.h>
#include <itensor/util/print_macro.h>
#include <vector>
#include <qcircuit.hpp>
#include <circuit_topology.hpp>
#include <quantum_gate.hpp>
#include <circuits.hpp>

using namespace qcircuit;

int main(int argc, char const* argv[]) {
    const auto topology = make_ibmq_topology();
    topology.exportDotString("test.dot");
    const size_t size = topology.numberOfBits();

    /* start with |00 ... 00> (53qubit) */
    std::vector<std::pair<std::complex<double>, std::complex<double>>>
        init_qbits(size, std::make_pair(1.0, 0.0));

    QCircuit circuit(topology, init_qbits);
    circuit.setCutoff(1e-5);

    /* Below is the demonstration of generating GHZ state */

    circuit.apply(H(6), X(11));    // apply Hadamard and X to gate (6,11)
    circuit.apply(H(10));          // apply Hadamard to gate 10
    circuit.apply(CNOT(10, 11));   // apply CNOT to gate (10, 11)
    circuit.apply(CNOT(6, 11));    // apply CNOT to gate (6, 11)
    circuit.apply(H(6), H(11));    // apply Hadamard to gate (6,11)
    circuit.apply(H(10));          // apply Hadamard to gate 10

    /* The result should be bell state (1/sqrt(2))(|000> + |111>) */

    /*
     * To show GHZ state is generated, calc the overlap between |0...000....0> and |0...111....0>,
     * where 000 and 111 are located on the qubit (6,10,11).
     */

    QCircuit circuit000(topology, init_qbits, circuit.site()); // |0...000....0>

    QCircuit circuit111(topology, init_qbits, circuit.site()); // to be |0...111....0> just below
    circuit111.apply(X(6), X(11)); // flip the qubit number (6,11)
    circuit111.apply(X(10));       // flip the qubit number 10

    std::vector<ITensor> op;
    op.reserve(size);
    for(size_t i = 0;i < size;i++) {
        op.push_back(circuit.generateTensorOp(Id(i)));
    }

    Print(overlap(circuit, op, circuit000)); //result should be -1/sqrt(2)
    Print(overlap(circuit, op, circuit111)); //result should be 1/sqrt(2)
    Print(overlap(circuit, op, circuit));    //result should be 1

    return 0;
}
```

## Python example
At the root directory, `pip install .` is available.
If you want to update existing one, `--no-cache-dir` option may be
required to avoid using old cache.

If your own ITensor is in a different place from `external/itensor`,
specify the location by editting `setup.py`.

### Python interface
```python
from qcircuit.core import *


topology = make_ibmq_topology()
circuit = QCircuit(topology)
circuit.cutoff = 1e-5

circuit.apply(H(0))
prob0 = circuit.probability_of_zero(0)
print("Probability to observe |0>: {:.3f}".format(prob0)) # should be 1/2

bit = circuit.observe_qubit(0)  # projection
print("Qubit 0 is observed as |{}>".format(bit))
```

### QASM interface
Currently under development.

```python
from qcircuit.qasm import *


data = """
OPENQASM 2.0;

include "qelib1.inc";

creg c[2];
qreg q[2];
x q[0];
measure q[0] -> c[0];
"""

engine = QASMInterpreter(data)
engine.execute()

print("{:02b}".format(engine._cregs.get_data("c")))
```
