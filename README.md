# IBM Q simulator 

An IBM Q simulator of 53-qubit quantum processor (rochester) by using tensor contraction (tensor network)

## Installation

Install [ITensor](https://itensor.org/) and change `LIBRARY_DIR` in Makefile to your ITensor directory.

## Example

```c++
#include "itensor/all.h"
#include "itensor/mps/siteset.h"
#include <stdio.h>
#include <fstream>
#include <cmath>
#include "ibmq.h"

using namespace itensor;
using namespace ibmq;

int main(int argc, char const* argv[]){

    std::array<std::pair<double, double>, IBMQPeps::NUM_BITS> init_param;

    //start with |0> (53qubit)
    for(auto& elem : init_param){
        elem = std::make_pair(1.0, 0.0);
    }

    IBMQPeps circuit(init_param);


    //Below is the demonstration of genearating GHZ state

    //the default cursor is located on qubit number 1 and 2
    //shift the cursor to qubit number (7,12)
    circuit.shift_to(3, {"Cutoff", 1E-5});
    circuit.shift_to(4, {"Cutoff", 1E-5});
    circuit.shift_to(5, {"Cutoff", 1E-5});
    circuit.shift_to(7, {"Cutoff", 1E-5});
    circuit.shift_to(12, {"Cutoff", 1E-5});

    //apply Hadamard and X to gate (7,12)
    circuit.apply(H(circuit.site(7))*X(circuit.site(12)));
    //shift the cursor to qubit number (11,12)
    circuit.shift_to(11, {"Cutoff", 1E-5});
    //apply Hadamard to gate 11
    circuit.apply(H(circuit.site(11))*Id(circuit.site(12)));
    //apply CNOT to gate (11, 12)
    circuit.apply(CNOT(circuit.site(11), circuit.site(12)));
    //shift the cursor to qubit number (7,12)
    circuit.shift_to(7, {"Cutoff", 1E-5});
    //apply CNOT to gate (7, 12)
    circuit.apply(CNOT(circuit.site(7), circuit.site(12)));
    //apply Hadamard to gate (7,12)
    circuit.apply(H(circuit.site(7))*H(circuit.site(12)));
    //shift the cursor to qubit number (11,12)
    circuit.shift_to(11, {"Cutoff", 1E-5});
    //apply Hadamard to gate 11
    circuit.apply(H(circuit.site(11))*Id(circuit.site(12)));
    //the result should be bell state (1/sqrt(2))(|000> + |111>)
    
    //shift the cursor to qubit number (7,12)
    circuit.shift_to(7, {"Cutoff", 1E-5});
    circuit.shift_to(5, {"Cutoff", 1E-5});
    circuit.shift_to(4, {"Cutoff", 1E-5});
    circuit.shift_to(3, {"Cutoff", 1E-5});
    circuit.shift_to(2, {"Cutoff", 1E-5});
    circuit.shift_to(1, {"Cutoff", 1E-5});


    //to show GHZ state is generated, calc the overlap between |0...000....0> and |0...111....0> where 000 and 111 are located on the qubit (7,11,12).
    
    //|0...000....0>
    IBMQPeps circuit000(init_param, circuit.site());

    //|0...111....0>
    IBMQPeps circuit111(init_param, circuit.site());
    //shift the cursor to qubit number (7,12)
    circuit111.shift_to(3, {"Cutoff", 1E-5});
    circuit111.shift_to(4, {"Cutoff", 1E-5});
    circuit111.shift_to(5, {"Cutoff", 1E-5});
    circuit111.shift_to(7, {"Cutoff", 1E-5});
    circuit111.shift_to(12, {"Cutoff", 1E-5});
    //flip the qubit number (7,12)
    circuit111.apply(X(circuit.site(7))*X(circuit.site(12)));
    //shift the cursor to qubit number (11,12)
    circuit111.shift_to(11, {"Cutoff", 1E-5});
    //flip the qubit number 11
    circuit111.apply(X(circuit.site(11))*Id(circuit.site(12)));
    
    myvector<ITensor> op(IBMQPeps::NUM_BITS);
    for(auto i : range1(IBMQPeps::NUM_BITS)){
        op[i] = Id(circuit.site(i));
    }

    Print(overlap(circuit, op, circuit000)); //result should be -1/sqrt(2)
    Print(overlap(circuit, op, circuit111)); //result should be 1/sqrt(2)
    Print(overlap(circuit, op, circuit)); //result should be 1

	return 0;
}
```


