import qiskit.qasm
import qiskit.qasm.node
from qcircuit.core import *
from .registers import QuantumRegisters, ClassicalRegisters


class QASMInterpreter:
    def __init__(self, data):
        self._data = data
        topology = make_ibmq_topology()
        self._engine = QCircuit(topology)
        self._max_qubit = topology.number_of_bits()
        self._qregs = QuantumRegisters(self._max_qubit)
        self._cregs = ClassicalRegisters()


    def execute(self):
        statements = qiskit.qasm.Qasm(data=self._data).parse().children
        for stat in statements:
            self._execute_statement(stat, {})


    def _execute_statement(self, stat, env):
        if type(stat) == qiskit.qasm.node.Format:
            pass # TODO
        elif type(stat) == qiskit.qasm.node.Qreg:
            name = stat.children[0].name
            size = stat.children[0].index
            self._qregs.add(name, size)
        elif type(stat) == qiskit.qasm.node.Creg:
            name = stat.children[0].name
            size = stat.children[0].index
            self._cregs.add(name, size)
        elif type(stat) == qiskit.qasm.node.Barrier:
            pass # TODO
        elif type(stat) == qiskit.qasm.node.Gate:
            self._define_gate(stat.children)
        elif type(stat) == qiskit.qasm.node.Measure:
            self._measure(stat.children)
        elif type(stat) == qiskit.qasm.node.Reset:
            self._reset(stat.children)
        elif type(stat) == qiskit.qasm.node.UniversalUnitary:
            self._call_universal_unitary(stat.children, env)
        elif type(stat) == qiskit.qasm.node.Cnot:
            self._call_cnot(stat.children)
        elif type(stat) == qiskit.qasm.node.CustomUnitary:
            self._call_custom_unitary(stat.children)
        elif type(stat) == qiskit.qasm.node.If:
            self._process_if_statement(stat.children)
        else:
            raise Exception('Unimplemented operation')


    def _define_gate(self, args):
        pass # TODO


    def _measure(self, args):
        if type(args[0]) == qiskit.qasm.node.IndexedId and type(args[1]) == qiskit.qasm.node.IndexedId:
            qubit_index = self._qregs.get_hardware_index(args[0].name, args[0].index)
            observed = self._engine.observe_qubit(qubit_index)
            self._cregs.set(args[1].name, args[1].index, observed)
        elif type(args[0]) == qiskit.qasm.node.Id and type(args[1]) == qiskit.qasm.node.Id:
            if self._qregs.get_size(args[0].name) != self._cregs.get_size(args[1].name):
                pass # TODO (throw Exception)

            size = self._qregs.get_size(args[0].name)
            for i in range(size):
                qubit_index = self._qregs.get_hardware_index(args[0].name, i)
                observed = self._engine.observe_qubit(qubit_index)
                self._cregs.set(args[1].name, i, observed)
        else:
            pass # TODO (throw Exception)


    def _reset(self, args):
        if type(args[0]) == qiskit.qasm.node.IndexedId:
            qubit_index = self._qregs.get_hardware_index(args[0].name, args[0].index)
            self._engine.reset_qubit(qreg_index)
        else:
            size = self._qregs.get_size(args[0].name)
            for i in range(size):
                qubit_index = self._qregs.get_hardware_index(args[0].name, i)
                self._engine.reset_qubit(qreg_index)


    def _call_universal_unitary(self, args, env):
        params = args[0].children
        qreg = args[1]
        if(args[0].size() != 3):
            pass # TODO (throw Exception)

        theta = params[0].sym([{}])
        phi = params[1].sym([{}])
        lamda = params[2].sym([{}]) # "lambda" is a reserved keyword

        if type(qreg) == qiskit.qasm.node.IndexedId:
            qubit_index = self._qregs.get_hardware_index(qreg.name, qreg.index)
            self._engine.apply(UniversalUnitary(qubit_index, theta, phi, lamda))
        elif type(qreg) == qiskit.qasm.node.Id:
            size = self._qregs.get_size(qreg.name)
            for i in range(size):
                qubit_index = self._qregs.get_hardware_index(qreg.name, i)
                self._engine.apply(UniversalUnitary(qubit_index, theta, phi, lamda))


    def _call_cnot(self, args):
        # CNOT must be able to receive any combination of indexed/unindexed registers

        if type(args[0]) == qiskit.qasm.node.IndexedId and type(args[1]) == qiskit.qasm.node.IndexedId:
            qubit_indices0 = [self._qregs.get_hardware_index(args[0].name, args[0].index)]
            qubit_indices1 = [self._qregs.get_hardware_index(args[1].name, args[1].index)]
        elif type(args[0]) == qiskit.qasm.node.IndexedId:
            size = self._qregs.get_size(args[1].name)
            qubit_indices0 = [self._qregs.get_hardware_index(args[0].name, args[0].index)] * size
            qubit_indices1 = [self._qregs.get_hardware_index(args[1].name, i) for i in range(size)]
        elif type(args[1]) == qiskit.qasm.node.IndexedId:
            size = self._qregs.get_size(args[0].name)
            qubit_indices0 = [self._qregs.get_hardware_index(args[0].name, i) for i in range(size)]
            qubit_indices1 = [self._qregs.get_hardware_index(args[1].name, args[1].index)] * size
        else:
            size = self._qregs.get_size(args[0].name)
            if size != self._qregs.get_size(args[1].name):
                pass # TODO (throw Exception)

            qubit_indices0 = [self._qregs.get_hardware_index(args[0].name, i) for i in range(size)]
            qubit_indices1 = [self._qregs.get_hardware_index(args[1].name, i) for i in range(size)]


        for i0, i1 in zip(qubit_indices0, qubit_indices1):
            self._engine.apply(CNOT(i0, i1))


    def _call_custom_unitary(self, args):
        pass # TODO


    def _process_if_statement(self, args):
        id = args[0].name
        value = args[1].value
        if id not in self._cregs:
            pass # TODO (throw Exception)

        if self._cregs.get_data(id) == value:
            self._execute_statement(args[2], {})
