#include <gtest/gtest.h>
#include <cmath>
#include <circuit_topology.hpp>


TEST(CIRCUIT_TOPOLOGY_TEST, CHECK_CONNECTED_GRAPH) {
    using namespace qcircuit;
    CircuitTopology topology(5);

    topology.generateLink(0, 1);
    topology.generateLink(0, 2);
    topology.generateLink(2, 3);
    topology.generateLink(3, 4);
    // 0 -+- 1
    //    +- 2 --- 3 --- 4

    EXPECT_TRUE(topology.isConnectedGraph());
}

TEST(CIRCUIT_TOPOLOGY_TEST, CHECK_NOT_CONNECTED_GRAPH) {
    using namespace qcircuit;
    CircuitTopology topology(5);

    topology.generateLink(0, 1);
    topology.generateLink(0, 2);
    topology.generateLink(3, 4);
    // 0 -+- 1
    //    +- 2     3 --- 4

    EXPECT_FALSE(topology.isConnectedGraph());
}
