# README for developers

## Dependencies
[ITensor](https://itensor.org/) v3 is required.
By default, CMakefile refers to `external/itensor` directory.

## Example

```c++
#include "itensor/all.h"
#include "itensor/util/print_macro.h"
#include <vector>
#include "qcircuit.hpp"
#include "circuit_topology.hpp"
#include "quantum_gate.hpp"
#include "circuits.hpp"

using namespace qcircuit;

int main(int argc, char const* argv[]){
    const auto topology = make_ibmq_topology();
    const size_t size = topology.numberOfBits();

    /* start with |00 ... 00> (53qubit) */
    std::vector<std::pair<std::complex<double>, std::complex<double>>>
        init_qbits(size, std::make_pair(1.0, 0.0));

    QCircuit circuit(topology, init_qbits);

    /* Below is the demonstration of generating GHZ state */

    circuit.apply(H(6), X(11), {"Cutoff", 1E-5});    // apply Hadamard and X to gate (6,11)
    circuit.apply(H(10), Id(11), {"Cutoff", 1E-5});  // apply Hadamard to gate 10
    circuit.apply(CNOT(10, 11), {"Cutoff", 1E-5});   // apply CNOT to gate (10, 11)
    circuit.apply(CNOT(6, 11), {"Cutoff", 1E-5});    // apply CNOT to gate (6, 11)
    circuit.apply(H(6), H(11), {"Cutoff", 1E-5});    // apply Hadamard to gate (6,11)
    circuit.apply(H(10), Id(11), {"Cutoff", 1E-5});  // apply Hadamard to gate 10

    /* The result should be bell state (1/sqrt(2))(|000> + |111>) */

    /*
     * To show GHZ state is generated, calc the overlap between |0...000....0> and |0...111....0>,
     * where 000 and 111 are located on the qubit (6,10,11).
     */

    QCircuit circuit000(topology, init_qbits, circuit.site()); // |0...000....0>

    QCircuit circuit111(topology, init_qbits, circuit.site()); // to be |0...111....0> just below
    circuit111.apply(X(6), X(11), {"Cutoff", 1E-5});   //flip the qubit number (6,11)
    circuit111.apply(X(10), Id(11), {"Cutoff", 1E-5}); //flip the qubit number 10

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


