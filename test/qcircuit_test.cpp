#include <gtest/gtest.h>
#include <cmath>
#include "itensor/util/print_macro.h"
#include "qcircuit.hpp"
#include "circuits.hpp"

TEST(CALCULATION_TEST, GHZ_STATE_TEST) {
    using namespace std;
    using namespace qcircuit;

    const auto topology = make_ibmq_topology();
    const size_t size = topology.numberOfBits();

    /* start with |00 ... 00> (53qubit) */
    vector<pair<complex<double>, complex<double>>>
        init_qbits(size, make_pair(1.0, 0.0));

    QCircuit circuit(topology, init_qbits);

    circuit.apply(H(6), X(11), {"Cutoff", 1E-5});
    circuit.apply(H(10), Id(11), {"Cutoff", 1E-5});
    circuit.apply(CNOT(10, 11), {"Cutoff", 1E-5});
    circuit.apply(CNOT(6, 11), {"Cutoff", 1E-5});
    circuit.apply(H(6), H(11), {"Cutoff", 1E-5});
    circuit.apply(H(10), Id(11), {"Cutoff", 1E-5});

    /* The result should be bell state (1/sqrt(2))(|000> + |111>) */

    /*
     * To show GHZ state is generated, calc the overlap between |0...000....0> and |0...111....0>,
     * where 000 and 111 are located on the qubit (6,10,11).
     */

    QCircuit circuit000(topology, init_qbits, circuit.site());  // |0...000....0>

    QCircuit circuit111(topology, init_qbits, circuit.site());  // to be |0...111....0> just below
    circuit111.apply(X(6), X(11), {"Cutoff", 1E-5});   // flip the qubit number (6,11)
    circuit111.apply(X(10), Id(11), {"Cutoff", 1E-5}); //flip the qubit number 10

    vector<ITensor> op;
    op.reserve(size);
    for(size_t i = 0;i < size;i++) {
        op.push_back(circuit.generateTensorOp(Id(i)));
    }

    EXPECT_NEAR(1/sqrt(2), abs(overlap(circuit, op, circuit000)), 1e-3);
    EXPECT_NEAR(1/sqrt(2), abs(overlap(circuit, op, circuit111)), 1e-3);
    EXPECT_NEAR(1.0 , abs(overlap(circuit, op, circuit)), 1e-3);
}


TEST(CALCULATION_TEST, LOOP_TEST) {
    /*  Almost the same as GHZ_STATE_TEST, but using a 8-site periodic chain. */
    using namespace std;
    using namespace qcircuit;

    const size_t size = 8;
    const auto topology = make_chain(size);

    vector<pair<complex<double>, complex<double>>>
        init_qbits(size, make_pair(1.0, 0.0));

    QCircuit circuit(topology, init_qbits);

    circuit.apply(H(0), X(1), {"Cutoff", 1E-5});
    circuit.apply(H(2), Id(1), {"Cutoff", 1E-5});
    circuit.apply(CNOT(2, 1), {"Cutoff", 1E-5});
    circuit.apply(Id(3), Id(4), {"Cutoff", 1E-5}); // make a detour
    circuit.apply(Id(5), Id(6), {"Cutoff", 1E-5}); //
    circuit.apply(Id(7), Id(0), {"Cutoff", 1E-5}); //
    circuit.apply(CNOT(0, 1), {"Cutoff", 1E-5});
    circuit.apply(H(0), H(1), {"Cutoff", 1E-5});
    circuit.apply(H(2), Id(1), {"Cutoff", 1E-5});

    QCircuit circuit000(topology, init_qbits, circuit.site());

    QCircuit circuit111(topology, init_qbits, circuit.site());
    circuit111.apply(X(0), X(1), {"Cutoff", 1E-5});
    circuit111.apply(X(2), Id(3), {"Cutoff", 1E-5});

    vector<ITensor> op;
    op.reserve(size);
    for(size_t i = 0;i < size;i++) {
        op.push_back(circuit.generateTensorOp(Id(i)));
    }

    EXPECT_NEAR(1/sqrt(2), abs(overlap(circuit, op, circuit000)), 1e-3);
    EXPECT_NEAR(1/sqrt(2), abs(overlap(circuit, op, circuit111)), 1e-3);
    EXPECT_NEAR(1.0 , abs(overlap(circuit, op, circuit)), 1e-3);

}
