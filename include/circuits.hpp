//Copyright (c) 2020 Jij Inc.


#pragma once

#include "circuit_topology.hpp"


namespace qcircuit {

    CircuitTopology make_ibmq_topology() {
        const size_t size = 53;
        CircuitTopology topology(size);

        topology.generateLink(0,1);
        topology.generateLink(1,2);
        topology.generateLink(2,3);
        topology.generateLink(3,4);

        topology.generateLink(0,5);
        topology.generateLink(4,6);
        topology.generateLink(5,7);
        topology.generateLink(6,11);

        topology.generateLink(7,8);
        topology.generateLink(8,9);
        topology.generateLink(9,10);
        topology.generateLink(10,11);

        topology.generateLink(7,12);
        topology.generateLink(11,13);
        topology.generateLink(12,14);
        topology.generateLink(13,15);
        topology.generateLink(14,16);
        topology.generateLink(15,18);

        topology.generateLink(9,17);

        topology.generateLink(16,19);
        topology.generateLink(18,20);
        topology.generateLink(19,21);
        topology.generateLink(20,22);
        topology.generateLink(21,23);
        topology.generateLink(22,27);

        topology.generateLink(17,25);

        topology.generateLink(23,24);
        topology.generateLink(24,25);
        topology.generateLink(25,26);
        topology.generateLink(26,27);

        topology.generateLink(23,28);
        topology.generateLink(27,29);
        topology.generateLink(28,30);
        topology.generateLink(29,34);

        topology.generateLink(30,31);
        topology.generateLink(31,32);
        topology.generateLink(32,33);
        topology.generateLink(33,34);

        topology.generateLink(30,35);
        topology.generateLink(34,36);
        topology.generateLink(35,37);
        topology.generateLink(36,38);
        topology.generateLink(37,39);
        topology.generateLink(38,41);

        topology.generateLink(32,40);

        topology.generateLink(39,42);
        topology.generateLink(41,43);
        topology.generateLink(42,44);
        topology.generateLink(43,45);
        topology.generateLink(44,46);
        topology.generateLink(45,50);

        topology.generateLink(40,48);

        topology.generateLink(46,47);
        topology.generateLink(47,48);
        topology.generateLink(48,49);
        topology.generateLink(49,50);

        topology.generateLink(46,51);
        topology.generateLink(50,52);

        return topology;
    }

    CircuitTopology make_chain(size_t size, bool periodic = true) {
        CircuitTopology topology(size);

        for(size_t i = 0;i < size-1;i++) {
            topology.generateLink(i, i+1);
        }

        if(periodic) {
            topology.generateLink(size-1, 0);
        }

        return topology;
    }

    CircuitTopology make_alltoall_topology(size_t size){
        CircuitTopology topology(size);
        for(size_t i = 0; i < size; i++){
            for(size_t j = i+1; j < size; j++){
                topology.generateLink(i, j);
            }
        }
        return topology;
    }
}
