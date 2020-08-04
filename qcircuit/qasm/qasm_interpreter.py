import qiskit.qasm
import qiskit.qasm.node
from qcircuit.core import *
from .registers import QuantumRegisters, ClassicalRegisters
from .exception import QASMError


class QASMInterpreter:
    def __init__(self, data):
        self._data = data
        topology = make_ibmq_topology()
        self._engine = QCircuit(topology)
        self._max_qubit = topology.number_of_bits()
        self._qregs = QuantumRegisters(self._max_qubit)
        self._cregs = ClassicalRegisters()
        self._custom_unitaries = {}
        self._global_env = {}


    def execute(self):
        statements = qiskit.qasm.Qasm(data=self._data).parse().children
        for stat in statements:
            self._execute_statement(stat, self._global_env)


    def _execute_statement(self, stat, env):
        if type(stat) == qiskit.qasm.node.Format:
            pass # TODO
        elif type(stat) == qiskit.qasm.node.Qreg:
            self._add_quantum_register(stat.children)
        elif type(stat) == qiskit.qasm.node.Creg:
            self._add_classical_register(stat.children)
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
            self._call_cnot(stat.children, env)
        elif type(stat) == qiskit.qasm.node.CustomUnitary:
            self._call_custom_unitary(stat.children, env)
        elif type(stat) == qiskit.qasm.node.If:
            self._process_if_statement(stat.children, env)
        else:
            raise QASMError(stat.type, 'Unimplemented operation')


    def _add_quantum_register(self, args):
        name = args[0].name
        size = args[0].index
        self._qregs.add(name, size)
        self._global_env[name] = name


    def _add_classical_register(self, args):
        name = args[0].name
        size = args[0].index
        self._cregs.add(name, size)
        self._global_env[name] = name


    def _define_gate(self, args):
        id = args[0].name

        if id not in self._custom_unitaries:
            self._custom_unitaries[id] = args
        else:
            raise QASMError('gate', 'Multiple definition of {}'.format(id))
            # Qiskit parser also raises the error, so this would never be raised


    def _measure(self, args):
        if type(args[0]) == qiskit.qasm.node.IndexedId and type(args[1]) == qiskit.qasm.node.IndexedId:
            qubit_index = self._qregs.get_hardware_index(args[0].name, args[0].index)
            observed = self._engine.observe_qubit(qubit_index)
            self._cregs.set(args[1].name, args[1].index, observed)
        elif type(args[0]) == qiskit.qasm.node.Id and type(args[1]) == qiskit.qasm.node.Id:
            if self._qregs.get_size(args[0].name) != self._cregs.get_size(args[1].name):
                message = 'Unmatching register size ({} and {})'.format(args[0].name, args[1].name)
                raise QASMError('measure', message)

            size = self._qregs.get_size(args[0].name)
            for i in range(size):
                qubit_index = self._qregs.get_hardware_index(args[0].name, i)
                observed = self._engine.observe_qubit(qubit_index)
                self._cregs.set(args[1].name, i, observed)
        else:
            raise QASMError('measure', 'Unsupported argument type')


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
        num_args = args[0].size()
        if num_args != 3:
            how = 'few' if num_args < 3 else 'many'
            raise QASMError('U', 'Too {} arguments (required 3, but {})'.format(how, num_args))

        theta = params[0].sym([env])
        phi = params[1].sym([env])
        lamda = params[2].sym([env]) # "lambda" is a reserved keyword

        reg_name = env[qreg.name] # resolve symbol in the current scope

        hardstart = self._qregs.get_hardware_index(reg_name, 0)

        if type(qreg) == qiskit.qasm.node.IndexedId:
            qubit_index = self._qregs.get_hardware_index(reg_name, qreg.index)
            self._engine.apply(UniversalUnitary(qubit_index, theta, phi, lamda))
        elif type(qreg) == qiskit.qasm.node.Id:
            size = self._qregs.get_size(reg_name)
            for i in range(size):
                qubit_index = self._qregs.get_hardware_index(reg_name, i)
                self._engine.apply(UniversalUnitary(qubit_index, theta, phi, lamda))


    def _call_cnot(self, args, env):
        # CNOT must be able to receive any combination of indexed/unindexed registers

        reg_name0 = env[args[0].name] # resolve symbol in the current scope
        reg_name1 = env[args[1].name] # resolve symbol in the current scope

        if type(args[0]) == qiskit.qasm.node.IndexedId and type(args[1]) == qiskit.qasm.node.IndexedId:
            qubit_indices0 = [self._qregs.get_hardware_index(reg_name0, args[0].index)]
            qubit_indices1 = [self._qregs.get_hardware_index(reg_name1, args[1].index)]
        elif type(args[0]) == qiskit.qasm.node.IndexedId:
            size = self._qregs.get_size(reg_name1)
            qubit_indices0 = [self._qregs.get_hardware_index(reg_name0, args[0].index)] * size
            qubit_indices1 = [self._qregs.get_hardware_index(reg_name1, i) for i in range(size)]
        elif type(args[1]) == qiskit.qasm.node.IndexedId:
            size = self._qregs.get_size(reg_name0)
            qubit_indices0 = [self._qregs.get_hardware_index(reg_name0, i) for i in range(size)]
            qubit_indices1 = [self._qregs.get_hardware_index(reg_name1, args[1].index)] * size
        else:
            size = self._qregs.get_size(reg_name0)
            if size != self._qregs.get_size(reg_name1):
                message = 'Unmatching register size ({} and {})'.format(reg_name0, reg_name1)
                raise QASMError('CX', message)

            qubit_indices0 = [self._qregs.get_hardware_index(reg_name0, i) for i in range(size)]
            qubit_indices1 = [self._qregs.get_hardware_index(reg_name1, i) for i in range(size)]


        for i0, i1 in zip(qubit_indices0, qubit_indices1):
            self._engine.apply(CNOT(i0, i1))


    def _call_custom_unitary(self, args, env):
        gate_info = self._custom_unitaries[args[0].name]

        reg_arg_pos = 2
        # The position of the register arguments.
        # This can be changed depending on whether the gate
        # receives c-number arguments or not.


        inner_env = {}
        if len(gate_info) > 3:
            # If the gate receives some c-number arguments

            for var, exp in zip(gate_info[1].children, args[1].children):
                inner_env[var.name] = qiskit.qasm.node.Real(exp.sym([env]))
        else:
            # If the gate receives no c-number arguments
            reg_arg_pos = 1


            # Gate receives at least one register argument.
            # (Because if not, what is the gate used for?)
        for reg_var, reg_name in zip(gate_info[reg_arg_pos].children, args[reg_arg_pos].children):
            inner_env[reg_var.name] = env[reg_name.name]

        for stat in gate_info[reg_arg_pos+1].children:
            self._execute_statement(stat, inner_env)


    def _process_if_statement(self, args, env):
        id = args[0].name
        value = args[1].value
        if id not in self._cregs:
            raise QASMError('if', 'Classical register {} not found'.format(id))

        if self._cregs.get_data(id) == value:
            self._execute_statement(args[2], env)
