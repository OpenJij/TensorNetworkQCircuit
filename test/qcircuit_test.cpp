#include <gtest/gtest.h>
#include <cmath>
#include <utility>
#include <itensor/util/print_macro.h>
#include <qcircuit.hpp>
#include <circuits.hpp>

TEST(QCIRCUIT_TEST, CHECK_INITIAL_CURSOR_POSITION) {
    using namespace qcircuit;
    CircuitTopology topology(4);

    topology.generateLink(0, 3);
    topology.generateLink(0, 2);
    topology.generateLink(0, 1);

    QCircuit circuit(topology);

    EXPECT_EQ((std::pair<size_t, size_t>(0, 1)), circuit.getCursor());
}

TEST(QCIRCUIT_TEST, INVALID_CURSOR_MOVE) {
    using namespace std;
    using namespace qcircuit;
    size_t size = 6;
    CircuitTopology topology(size);
    topology.generateLink(0, 1); //
    topology.generateLink(0, 2); //
    topology.generateLink(0, 3); //
    topology.generateLink(0, 4); //
    topology.generateLink(0, 5); // Starfish-like structure

    QCircuit circuit(topology);

    circuit.setCutoff(1e-5);

    EXPECT_THROW( {
        circuit.apply(X(2), X(3)); // Not connected sites
    }, QCircuitException);
}

TEST(CALCULATION_TEST, GHZ_STATE_TEST) {
    using namespace std;
    using namespace qcircuit;

    const auto topology = make_ibmq_topology();
    const size_t size = topology.numberOfBits();

    /* start with |00 ... 00> (53qubit) */
    vector<pair<complex<double>, complex<double>>>
        init_qbits(size, make_pair(1.0, 0.0));

    QCircuit circuit(topology, init_qbits);
    circuit.setCutoff(1e-5);

    circuit.apply(H(6), X(11));
    circuit.apply(H(10));
    circuit.apply(CNOT(10, 11));
    circuit.apply(CNOT(6, 11));
    circuit.apply(H(6), H(11));
    circuit.apply(H(10));

    /* The result should be bell state (1/sqrt(2))(|000> + |111>) */

    /*
     * To show GHZ state is generated, calc the overlap between |0...000....0> and |0...111....0>,
     * where 000 and 111 are located on the qubit (6,10,11).
     */

    QCircuit circuit000(topology, init_qbits, circuit.site());  // |0...000....0>
    circuit000.setCutoff(1e-5);

    QCircuit circuit111(topology, init_qbits, circuit.site());  // to be |0...111....0> just below
    circuit111.setCutoff(1e-5);
    circuit111.apply(X(6), X(11)); // flip the qubit number (6,11)
    circuit111.apply(X(10));       // flip the qubit number 10

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
    circuit.setCutoff(1e-5);

    circuit.apply(H(0), X(1));
    circuit.apply(H(2), Id(1));
    circuit.apply(CNOT(2, 1));
    circuit.moveCursorAlong({3, 4, 5, 6, 7, 0}); // make a detour
    circuit.apply(CNOT(0, 1));
    circuit.apply(H(0), H(1));
    circuit.apply(H(2), Id(1));

    std::cout << "+++++++++++++ end cursor move +++++++++++++" << std::endl;

    QCircuit circuit000(topology, init_qbits, circuit.site());
    circuit000.setCutoff(1e-5);

    QCircuit circuit111(topology, init_qbits, circuit.site());
    circuit111.setCutoff(1e-5);
    circuit111.apply(X(0), X(1));
    circuit111.apply(X(2), Id(3));

    vector<ITensor> op;
    op.reserve(size);
    for(size_t i = 0;i < size;i++) {
        op.push_back(circuit.generateTensorOp(Id(i)));
    }

    EXPECT_NEAR(1/sqrt(2), abs(overlap(circuit, op, circuit000)), 1e-3);
    EXPECT_NEAR(1/sqrt(2), abs(overlap(circuit, op, circuit111)), 1e-3);
    EXPECT_NEAR(1.0 , abs(overlap(circuit, op, circuit)), 1e-3);

}

TEST(CALCULATION_TEST, SWAP_TEST) {
    using namespace std;
    using namespace qcircuit;

    const size_t size = 8;
    const auto topology = make_chain(size);

    vector<pair<complex<double>, complex<double>>>
        init_qbits(size, make_pair(1.0, 0.0));

    QCircuit circuit(topology, init_qbits);
    circuit.setCutoff(1e-5);
    circuit.apply(Id(0), X(1));
    circuit.apply(Swap(0, 1));

    QCircuit circuit10(topology, init_qbits, circuit.site());
    circuit10.setCutoff(1e-5);
    circuit10.apply(X(0), Id(1));

    vector<ITensor> op;
    op.reserve(size);
    for(size_t i = 0;i < size;i++) {
        op.push_back(circuit.generateTensorOp(Id(i)));
    }

    EXPECT_NEAR(1, abs(overlap(circuit, op, circuit10)), 1e-3);
}

TEST(CALCULATION_TEST, OBSERVATION_TEST) {
    using namespace std;
    using namespace qcircuit;

    const size_t size = 8;
    const auto topology = make_chain(size);

    vector<pair<complex<double>, complex<double>>>
        init_qbits(size, make_pair(1.0, 0.0));

    QCircuit circuit(topology, init_qbits);
    circuit.apply(H(0), Id(1)); // create (|0>+|1>)/sqrt(2)
    auto zero_weight = circuit.probabilityOfZero(0);
    EXPECT_NEAR(0.5, zero_weight, 1e-3); // probability to observe 0 should be 1/2.

    auto bit = circuit.observeQubit(0);
    std::cout << "qubit is observed as " << bit << std::endl;

    vector<ITensor> op;
    op.reserve(size);
    for(size_t i = 0;i < size;i++) {
        op.push_back(circuit.generateTensorOp(Id(i)));
    }
    if(bit == 0) {
        std::cout << "checking overlap with |0>" << std::endl;
        QCircuit circuit0(topology, init_qbits, circuit.site());
        EXPECT_NEAR(1, abs(overlap(circuit, op, circuit0)), 1e-3);
    } else {
        std::cout << "checking overlap with |1>" << std::endl;
        QCircuit circuit1(topology, init_qbits, circuit.site());
        circuit1.apply(X(0), Id(1));
        EXPECT_NEAR(1, abs(overlap(circuit, op, circuit1)), 1e-3);
    }
}

TEST(CALCULATION_TEST, CIRCUIT_WITH_MULTIPLE_LINKS) {
    using namespace std;
    using namespace qcircuit;
    size_t size = 6;
    CircuitTopology topology(size);
    topology.generateLink(0, 1);
    topology.generateLink(0, 2);
    topology.generateLink(0, 3);
    topology.generateLink(0, 4);
    topology.generateLink(0, 5);

    vector<pair<complex<double>, complex<double>>>
        init_qbits(size, make_pair(1.0, 0.0));

    // generate maximal entangled state
    QCircuit circuit(topology, init_qbits);
    circuit.setCutoff(1e-5);
    circuit.apply(H(0));
    cout << "cursor: " << circuit.getCursor().first << " " << circuit.getCursor().second << endl;
    circuit.apply(CNOT(0,1));
    cout << "cursor: " << circuit.getCursor().first << " " << circuit.getCursor().second << endl;
    circuit.apply(CNOT(0,2));
    cout << "cursor: " << circuit.getCursor().first << " " << circuit.getCursor().second << endl;
    circuit.apply(CNOT(0,3));
    cout << "cursor: " << circuit.getCursor().first << " " << circuit.getCursor().second << endl;
    circuit.apply(CNOT(0,4));
    cout << "cursor: " << circuit.getCursor().first << " " << circuit.getCursor().second << endl;
    circuit.apply(CNOT(0,5));
    cout << "cursor: " << circuit.getCursor().first << " " << circuit.getCursor().second << endl;

    // |000000>
    QCircuit circuit000000(topology, init_qbits, circuit.site());
    circuit000000.setCutoff(1e-5);

    // |111111>
    QCircuit circuit111111(topology, init_qbits, circuit.site());
    circuit111111.setCutoff(1e-5);
    circuit111111.apply(X(0));
    circuit111111.apply(X(1));
    circuit111111.apply(X(2));
    circuit111111.apply(X(3));
    circuit111111.apply(X(4));
    circuit111111.apply(X(5));

    // the following code causes segfault
    // circuit111111.apply(X(0), X(1));
    // circuit111111.apply(X(2), X(3));
    // circuit111111.apply(X(4), X(5));

    //identity operator
    vector<ITensor> op;
    op.reserve(size);
    for(size_t i = 0;i < size;i++) {
        op.push_back(circuit.generateTensorOp(Id(i)));
    }

    EXPECT_NEAR(1/sqrt(2), abs(overlap(circuit, op, circuit000000)), 1e-3);
    EXPECT_NEAR(1/sqrt(2), abs(overlap(circuit, op, circuit111111)), 1e-3);
    EXPECT_NEAR(1.0 , abs(overlap(circuit, op, circuit)), 1e-3);
}

TEST(CALCULATION_TEST, CIRCUIT_WITH_ALL_TO_ALL_CONNECTIVITY) {
    using namespace std;
    using namespace qcircuit;
    size_t size = 6;
    CircuitTopology topology = make_alltoall_topology(size);

    vector<pair<complex<double>, complex<double>>>
        init_qbits(size, make_pair(1.0, 0.0));

    // generate maximal entangled state
    QCircuit circuit(topology, init_qbits);
    circuit.setCutoff(1e-5);
    circuit.apply(H(0));
    cout << "cursor: " << circuit.getCursor().first << " " << circuit.getCursor().second << endl;
    circuit.apply(CNOT(0,1));
    cout << "cursor: " << circuit.getCursor().first << " " << circuit.getCursor().second << endl;
    circuit.apply(CNOT(0,2));
    cout << "cursor: " << circuit.getCursor().first << " " << circuit.getCursor().second << endl;
    circuit.apply(CNOT(0,3));
    cout << "cursor: " << circuit.getCursor().first << " " << circuit.getCursor().second << endl;
    circuit.apply(CNOT(0,4));
    cout << "cursor: " << circuit.getCursor().first << " " << circuit.getCursor().second << endl;
    circuit.apply(CNOT(0,5));
    cout << "cursor: " << circuit.getCursor().first << " " << circuit.getCursor().second << endl;

    // |000000>
    QCircuit circuit000000(topology, init_qbits, circuit.site());
    circuit000000.setCutoff(1e-5);

    // |111111>
    QCircuit circuit111111(topology, init_qbits, circuit.site());
    circuit111111.setCutoff(1e-5);
    cout << "cursor: " << circuit111111.getCursor().first << " " << circuit111111.getCursor().second << endl;
    circuit111111.apply(X(0));                                                   
    cout << "cursor: " << circuit111111.getCursor().first << " " << circuit111111.getCursor().second << endl;
    circuit111111.apply(X(1));                                                   
    cout << "cursor: " << circuit111111.getCursor().first << " " << circuit111111.getCursor().second << endl;
    circuit111111.apply(X(2));                                                   
    cout << "cursor: " << circuit111111.getCursor().first << " " << circuit111111.getCursor().second << endl;
    circuit111111.apply(X(3));                                                   
    cout << "cursor: " << circuit111111.getCursor().first << " " << circuit111111.getCursor().second << endl;
    circuit111111.apply(X(4));                                                   
    cout << "cursor: " << circuit111111.getCursor().first << " " << circuit111111.getCursor().second << endl;
    circuit111111.apply(X(5));                                                   
    cout << "cursor: " << circuit111111.getCursor().first << " " << circuit111111.getCursor().second << endl;

    //identity operator
    vector<ITensor> op;
    op.reserve(size);
    for(size_t i = 0;i < size;i++) {
        op.push_back(circuit.generateTensorOp(Id(i)));
    }

    EXPECT_NEAR(1/sqrt(2), abs(overlap(circuit, op, circuit000000)), 1e-3);
    EXPECT_NEAR(1/sqrt(2), abs(overlap(circuit, op, circuit111111)), 1e-3);
    EXPECT_NEAR(1.0 , abs(overlap(circuit, op, circuit)), 1e-3);
}


TEST(CALCULATION_TEST, MULTIPLE_MEASUREMENT_TEST) {
    using namespace std;
    using namespace qcircuit;

    const size_t size = 8;
    const auto topology = make_chain(size);

    QCircuit circuit(topology);
    circuit.apply(Id(2));
    circuit.apply(Id(3));


    double prob0 = circuit.probabilityOf(2, 0);
    double prob1 = circuit.probabilityOf(2, 1);
    std::cout << "weight of |0> at 2: " << prob0 << std::endl;
    std::cout << "weight of |1> at 2: " << prob1 << std::endl;
    EXPECT_NEAR(1.0, prob0 + prob1, 1e-3);
    circuit.observeQubit(2);

    prob0 = circuit.probabilityOf(3, 0);
    prob1 = circuit.probabilityOf(3, 1);
    std::cout << "weight of |0> at 3: " << prob0 << std::endl;
    std::cout << "weight of |1> at 3: " << prob1 << std::endl;
    EXPECT_NEAR(1.0, prob0 + prob1, 1e-3);
    circuit.observeQubit(3);
}
