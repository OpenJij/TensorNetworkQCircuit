from qcircuit import *

def main():
    topology = make_ibmq_topology()
    circuit = QCircuit(topology)
    circuit.cutoff = 1e-5
    print("Cutoff is {:e}".format(circuit.cutoff))

    path = [2, 3, 4, 6, 11, 10, 9, 8, 7]
    circuit.moveCursorAlong(path)
    print("detour along : ", end="")
    print(*path)

    print("Each qubit are initialized with |0>")

    prob0 = circuit.probabilityOfZero(0)
    print("Probability to observe |0>: {:f}".format(prob0))

    circuit.apply(X(0));
    print("Applying X to qubit 0")

    prob0 = circuit.probabilityOfZero(0)
    print("Probability to observe |0>: {:f}".format(prob0))

    circuit.apply(H(0));
    print("Applying H to qubit 0")

    prob0 = circuit.probabilityOfZero(0)
    print("Probability to observe |0>: {:f}".format(prob0))

    bit = circuit.observeQubit(0)
    print("Qubit 0 is observed as |{}>".format(bit))

    prob0 = circuit.probabilityOfZero(0)
    print("Probability to observe |0>: {:f}".format(prob0))


if __name__ == "__main__":
    main()
