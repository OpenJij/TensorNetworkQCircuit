import qiskit.qasm
import qiskit.qasm.node
from qcircuit.core import *

class QASMInterpreter:
    def __init__(self, data):
        self._data = data
        topology = make_ibmq_topology()
        self._engine = QCircuit(topology)
        self._max_qubit = topology.number_of_bits()
        self._creg = [0] * self._max_qubit
        self._qreg_top = 0
        self._qreg_info = {}
        self._creg_top = 0
        self._creg_info = {}


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
            self._add_qreg(name, size)
        elif type(stat) == qiskit.qasm.node.Creg:
            name = stat.children[0].name
            size = stat.children[0].index
            self._add_creg(name, size)
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
        else:
            raise Exception('Unimplemented operation')


    def _add_qreg(self, name, size):
        if self._qreg_top + size > self._max_qubit:
            pass # TODO (throw exception)

        self._qreg_info[name] = RegisterInfo(self._qreg_top, size)
        self._qreg_top += size


    def _get_qreg_index(self, name, index):
        qreg = self._qreg_info[name]
        if index <= qreg.size:
            pass # TODO (throw exception)

        return qreg.start + index


    def _get_qreg_size(self, name):
        return self._qreg_info[name].size


    def _add_creg(self, name, size):
        self._creg_info[name] = RegisterInfo(self._creg_top, size)
        self._creg_top += size


    def _get_creg_index(self, name, index):
        creg = self._creg_info[name]
        if index <= creg.size:
            pass # TODO (throw exception)

        return creg.start + index


    def _get_creg_size(self, name):
        return self._creg_info[name].size


    def _define_gate(self, args):
        pass # TODO


    def _measure(self, args):
        if type(args[0]) == qiskit.qasm.node.IndexedId and type(args[1]) == qiskit.qasm.node.IndexedId:
            qreg_index = self._get_qreg_index(args[0].name, args[0].index)
            creg_index = self._get_creg_index(args[1].name, args[1].index)

            self._creg[creg_index] = self._engine.observe_qubit(qreg_index)
        elif type(args[0]) == qiskit.qasm.node.Id and type(args[1]) == qiskit.qasm.node.Id:
            if self._get_qreg_size(args[0].name) != self._get_creg_size(args[1].name):
                pass # TODO (throw Exception)

            size = self._get_qreg_size(args[1].name)
            for i in range(size):
                qreg_index = self._get_qreg_index(args[0].name, args[0].index + i)
                creg_index = self._get_creg_index(args[1].name, args[1].index + i)

                self._creg[creg_index] = self._engine.observe_qubit(qreg_index)
        else:
            pass # TODO (throw Exception)


    def _reset(self, args):
        pass # TODO


    def _call_universal_unitary(self, args, env):
        # TODO:
        # Universal unitary must be able to receive unindexed register

        params = args[0].children
        if(args[0].size() != 3):
            pass # TODO (throw Exception)

        theta = params[0].sym([{}])
        phi = params[1].sym([{}])
        lamda = params[2].sym([{}]) # "lambda" is a reserved keyword

        qreg_index = self._get_qreg_index(args[1].name, args[1].index)

        self._engine.apply(UniversalUnitary(qreg_index, theta, phi, lamda));


    def _call_cnot(self, args):
        if type(args[0]) != qiskit.qasm.node.IndexedId or type(args[1]) != qiskit.qasm.node.IndexedId:
            pass # TODO (throw Exception)

        qreg_index0 = _get_qreg_index(args[0].name, args[0].index)
        qreg_index1 = _get_qreg_index(args[1].name, args[1].index)

        self._engine.apply(CNOT(qreg_index0, qreg_index1))


    def _call_custom_unitary(self, args):
        pass # TODO


class RegisterInfo:
    def __init__(self, start, size):
        self.start = start
        self.size = size

    def __repr__(self):
        return "RegisterInfo(start={}, size={})".format(self.start, self.size)
