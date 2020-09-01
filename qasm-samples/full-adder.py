# Quantum 2-bit full adder

from qcircuit.qasm import *

data ="""
OPENQASM 2.0;

include "qelib1.inc";

qreg q[4];
creg c[4];

h q[2]; // input
h q[3]; // input

ccx q[2], q[3], q[1]; // carry
cx q[2], q[0]; // sum
cx q[3], q[0]; // sum

measure q -> c;
"""

result = {}

for i in range(100):
    print("iteration #{}".format(i))
    engine = QASMInterpreter(data)
    engine.execute()

    bits = "{:04b}".format(engine._cregs.get_data("c"))

    result[bits] = result.get(bits, 0) + 1 # count appearance

print("Upper (lefter) two digits are input and lower (righter) two digit are result:")
print(result)
