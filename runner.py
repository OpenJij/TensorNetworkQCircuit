from qcircuit import *

def main():
    topology = make_ibmq_topology()
    circuit = QCircuit(topology)
    print("Each qubit are initialized with |0>")

    prob0 = circuit.probabilityOfZero(0)
    print("Probability to observe |0>: {:.3f}".format(prob0))

    circuit.apply(X(0), Id(1))
    print("Applying X to qubit 0")

    prob0 = circuit.probabilityOfZero(0)
    print("Probability to observe |0>: {:.3f}".format(prob0))

    circuit.apply(H(0), Id(1))
    print("Applying H to qubit 0")

    prob0 = circuit.probabilityOfZero(0)
    print("Probability to observe |0>: {:.3f}".format(prob0))

    bit = circuit.observeQubit(0)
    print("Qubit 0 is observed as |{}>".format(bit))

    prob0 = circuit.probabilityOfZero(0)
    print("Probability to observe |0>: {:.3f}".format(prob0))


if __name__ == "__main__":
    main()
