//Copyright (c) 2019 Jij Inc.

#include "itensor/all.h"
#include "itensor/util/print_macro.h"
#include <complex>
#include <vector>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include "qcircuit.hpp"
#include "circuit_topology.hpp"
#include "quantum_gate.hpp"

using namespace qcircuit;

CircuitTopology make_ibmq_topology();

int main(int argc, char const* argv[]){
    const auto topology = make_ibmq_topology();
    const size_t size = topology.numberOfBits();

    //start with |0> (53qubit)
    std::vector<std::pair<std::complex<double>, std::complex<double>>>
        init_qbits(size, std::make_pair(1.0, 0.0));

    QCircuit circuit(topology, init_qbits);

    //Below is the demonstration of generating GHZ state

    //the default cursor is located on qubit number 0 and 1

    //apply Hadamard and X to gate (6,11)
    circuit.apply(H(6), X(11), {"Cutoff", 1E-5});

    //apply Hadamard to gate 10
    circuit.apply(H(10), Id(11), {"Cutoff", 1E-5});

    //apply CNOT to gate (10, 11)
    circuit.apply(CNOT(10, 11), {"Cutoff", 1E-5});

    //apply CNOT to gate (6, 11)
    circuit.apply(CNOT(6, 11), {"Cutoff", 1E-5});

    //apply Hadamard to gate (6,11)
    circuit.apply(H(6), H(11), {"Cutoff", 1E-5});

    //apply Hadamard to gate 10
    circuit.apply(H(10), Id(11), {"Cutoff", 1E-5});

    //the result should be bell state (1/sqrt(2))(|000> + |111>)

    //to show GHZ state is generated, calc the overlap between |0...000....0> and |0...111....0> where 000 and 111 are located on the qubit (7,11,12).

    //|0...000....0>
    QCircuit circuit000(topology, init_qbits, circuit.site());

    //|0...111....0>
    QCircuit circuit111(topology, init_qbits, circuit.site());
    //flip the qubit number (6,11)
    circuit111.apply(X(6), X(11));

    //flip the qubit number 10
    circuit111.apply(X(10), Id(11));

    std::vector<ITensor> op;
    op.reserve(size);
    for(size_t i = 0;i < size;i++) {
        op.push_back(circuit.generateTensorOp(Id(i)));
    }

    Print(overlap(circuit, op, circuit000)); //result should be -1/sqrt(2)
    Print(overlap(circuit, op, circuit111)); //result should be 1/sqrt(2)
    Print(overlap(circuit, op, circuit)); //result should be 1

    return 0;
}



CircuitTopology make_ibmq_topology() {
    const size_t size = 53;
    CircuitTopology topology(size);

    topology.generate_link(0,1);
    topology.generate_link(1,2);
    topology.generate_link(2,3);
    topology.generate_link(3,4);

    topology.generate_link(0,5);
    topology.generate_link(4,6);
    topology.generate_link(5,7);
    topology.generate_link(6,11);

    topology.generate_link(7,8);
    topology.generate_link(8,9);
    topology.generate_link(9,10);
    topology.generate_link(10,11);

    topology.generate_link(7,12);
    topology.generate_link(11,13);
    topology.generate_link(12,14);
    topology.generate_link(13,15);
    topology.generate_link(14,16);
    topology.generate_link(15,18);

    topology.generate_link(9,17);

    topology.generate_link(16,19);
    topology.generate_link(18,20);
    topology.generate_link(19,21);
    topology.generate_link(20,22);
    topology.generate_link(21,23);
    topology.generate_link(22,27);

    topology.generate_link(17,25);

    topology.generate_link(23,24);
    topology.generate_link(24,25);
    topology.generate_link(25,26);
    topology.generate_link(26,27);

    topology.generate_link(23,28);
    topology.generate_link(27,29);
    topology.generate_link(28,30);
    topology.generate_link(29,34);

    topology.generate_link(30,31);
    topology.generate_link(31,32);
    topology.generate_link(32,33);
    topology.generate_link(33,34);

    topology.generate_link(30,35);
    topology.generate_link(34,36);
    topology.generate_link(35,37);
    topology.generate_link(36,38);
    topology.generate_link(37,39);
    topology.generate_link(38,41);

    topology.generate_link(32,40);

    topology.generate_link(39,42);
    topology.generate_link(41,43);
    topology.generate_link(42,44);
    topology.generate_link(43,45);
    topology.generate_link(44,46);
    topology.generate_link(45,50);

    topology.generate_link(40,48);

    topology.generate_link(46,47);
    topology.generate_link(47,48);
    topology.generate_link(48,49);
    topology.generate_link(49,50);

    topology.generate_link(46,51);
    topology.generate_link(50,52);

    return topology;
}
