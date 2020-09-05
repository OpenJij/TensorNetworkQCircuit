import qiskit.qasm
import qiskit.qasm.node
from qcircuit.core import *
from .registers import QuantumRegisters, ClassicalRegisters
from .exception import QASMError


class QASMInterpreter:
    """
    The QASMInterpreter engine, supporing OpenQASM 2.0.

    Attributes:
        _data (string): The QASM code.
        _engine (QCircuit): The quantum circuit object.
        _max_qubit (int): The number of hardware qubits.
        _qregs (QuautumRegisters): The quantum register information.
        _cregs (ClassicalRegisters): The classical register information.
        _custom_unitaries (dict[str, list[qiskit.qasm.node.Node]]): The dictionary which maps
            a function name to its definition.
    """

    def __init__(self, data, topology = make_chain(50)):
        """
        Constructs QASM interpreter.

        Even if invalid QASM code is given, this constructor raises NO errors.

        Args:
            data (str): The QASM code.
            topology (CircuitTopology): The qauntum circuit topology.
                Default is a one dimensional chain consisting of 50 qubits.
        """

        self._data = data
        self._engine = QCircuit(topology)
        self._max_qubit = topology.number_of_bits()
        self._qregs = QuantumRegisters(self._max_qubit)
        self._cregs = ClassicalRegisters()
        self._custom_unitaries = {}

    def execute(self):
        """
        Executes registered QASM code.

        Raises:
            QASMError: The registered QASM code is invalid.
        """

        try:
            statements = qiskit.qasm.Qasm(data=self._data).parse().children
        except qiskit.qasm.exceptions.QasmError as ex:
            raise QASMError('Parsing QASM', str(ex))

        for stat in statements:
            self._execute_statement(stat, {})

    def _execute_statement(self, stat, env):
        """
        Executes a single statement inside given variable environment.

        Args:
            stat (qiskit.qasm.node.Node): The node to be evaluated.
            env (dict[str, qiskit.qasm.node.Node): The variable environment in which
                variable are evaluated to classical values or QASM register IDs.

        Raises:
            QASMError: The registered QASM code is invalid or some unimplemented operation is called.
        """

        if type(stat) == qiskit.qasm.node.Format:
            pass  # TODO
        elif type(stat) == qiskit.qasm.node.Qreg:
            self._add_quantum_register(stat.children)
        elif type(stat) == qiskit.qasm.node.Creg:
            self._add_classical_register(stat.children)
        elif type(stat) == qiskit.qasm.node.Barrier:
            pass  # TODO
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
        """
        Registers a quantum register.

        Args:
            args (list[qiskit.qasm.node.Node]): The quantum register info to be registered.

        Raises:
            QASMError: The registration is failed.
        """
        name = args[0].name
        size = args[0].index
        self._qregs.add(name, size)

    def _add_classical_register(self, args):
        """
        Registers a classical register.

        Args:
            args (list[qiskit.qasm.node.Node]): The classical register info to be registered.

        Raises:
            QASMError: The registration is failed.
        """

        name = args[0].name
        size = args[0].index
        self._cregs.add(name, size)

    def _define_gate(self, args):
        """
        Defines a quantum gate.

        Args:
            args (list[qiskit.qasm.node.Node]): The quantum gate info to be registered.

        Raises:
            QASMError: The gate name is already used
                (This would never be raised because the same check is performed at parsing phase).
        """
        id = args[0].name

        if id not in self._custom_unitaries:
            self._custom_unitaries[id] = args
        else:
            raise QASMError('gate', 'Multiple definition of {}'.format(id))
            # This error would never be raised because
            # Qiskit parser also raises the error
            # wtih multiple definition

    def _measure(self, args):
        """
        Measures qubit(s) and stores the result into a classical register.

        Args:
            args (list[qiskit.qasm.node.Node]): The quantum register and
                classical register info to be used in the measuremnt.

        Raises:
            QASMError: Specifying invalid combination of quantum and classical registers.
        """
        if (type(args[0]) == qiskit.qasm.node.IndexedId and
                type(args[1]) == qiskit.qasm.node.IndexedId):
            qubit_index = self._qregs.get_hardware_index(args[0].name, args[0].index)
            observed = self._engine.observe_qubit(qubit_index)
            self._cregs.set(args[1].name, args[1].index, observed)
        elif (type(args[0]) == qiskit.qasm.node.Id and
                type(args[1]) == qiskit.qasm.node.Id):
            if (self._qregs.get_size(args[0].name) != self._cregs.get_size(args[1].name)):
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
        """
        Reset qubit(s) to |0>.

        Args:
            args (list[qiskit.qasm.node.Node]): The quantum register info to be reset

        Note:
           The OpenQASM 2.0 specification requires the reset operation to generates a mixed state.
           However, the spec can't be implemented due to an engine limitation.
           Thus the current version of this method just makes projection onto |0>.
        """

        if type(args[0]) == qiskit.qasm.node.IndexedId:
            qubit_index = self._qregs.get_hardware_index(args[0].name, args[0].index)
            self._engine.reset_qubit(args[0].name)
        else:
            size = self._qregs.get_size(args[0].name)
            for i in range(size):
                qubit_index = self._qregs.get_hardware_index(args[0].name, i)
                self._engine.reset_qubit(qubit_index)

    def _call_universal_unitary(self, args, env):
        """
        Calls a universal unitary gate.

        Args:
            args (list[qiskit.qasm.node.Node]): The universal unitary gate info to be called.
            env (dict[str, qiskit.qasm.node.Node): The variable environment in which
                variable are evaluated to classical values or QASM register IDs.
        """

        params = args[0].children
        num_args = args[0].size()
        if num_args != 3:
            how = 'few' if num_args < 3 else 'many'
            raise QASMError('U', 'Too {} arguments (required 3, but {})'.format(how, num_args))

        theta = params[0].sym([env])
        phi = params[1].sym([env])
        lamda = params[2].sym([env])  # "lambda" is a reserved keyword

        if not env:
            qreg = args[1]
        else:
            qreg = env[args[1].name]  # resolve symbol in the current scope

        if type(qreg) == qiskit.qasm.node.IndexedId:
            qubit_index = self._qregs.get_hardware_index(qreg.name, qreg.index)
            self._engine.apply(UniversalUnitary(qubit_index, theta, phi, lamda))
        elif type(qreg) == qiskit.qasm.node.Id:
            size = self._qregs.get_size(qreg.name)
            for i in range(size):
                qubit_index = self._qregs.get_hardware_index(qreg.name, i)
                self._engine.apply(UniversalUnitary(qubit_index, theta, phi, lamda))

    def _call_cnot(self, args, env):
        """
        Calls a Controlled-Not gate.

        Args:
            args (list[qiskit.qasm.node.Node]): The CNOT gate info to be called.
            env (dict[str, qiskit.qasm.node.Node): The variable environment in which
                variable are evaluated to classical values or QASM register IDs.

        Raises:
            QASMError: Specifying different size of quantum registers.
        """

        # CNOT must be able to receive any combination of
        # indexed/unindexed registers

        if not env:
            qreg0 = args[0]
            qreg1 = args[1]
        else:
            qreg0 = env[args[0].name]  # resolve symbol in the current scope
            qreg1 = env[args[1].name]  # resolve symbol in the current scope

        if type(qreg0) == qiskit.qasm.node.IndexedId and type(qreg1) == qiskit.qasm.node.IndexedId:
            qubit_v_indices0 = [self._qregs.get_virtual_index(qreg0.name, qreg0.index)]
            qubit_v_indices1 = [self._qregs.get_virtual_index(qreg1.name, qreg1.index)]
        elif type(qreg0) == qiskit.qasm.node.IndexedId:
            size = self._qregs.get_size(qreg1.name)
            qubit_v_indices0 = [self._qregs.get_virtual_index(qreg0.name, qreg0.index)] * size
            qubit_v_indices1 = [self._qregs.get_virtual_index(qreg1.name, i) for i in range(size)]
        elif type(qreg1) == qiskit.qasm.node.IndexedId:
            size = self._qregs.get_size(qreg0.name)
            qubit_v_indices0 = [self._qregs.get_virtual_index(qreg0.name, i) for i in range(size)]
            qubit_v_indices1 = [self._qregs.get_virtual_index(qreg1.name, qreg1.index)] * size
        else:
            size = self._qregs.get_size(qreg0.name)
            if size != self._qregs.get_size(qreg1.name):
                message = 'Unmatching register size ({} and {})'.format(qreg0.name, qreg1.name)
                raise QASMError('CX', message)

            qubit_v_indices0 = [self._qregs.get_virtual_index(qreg0.name, i) for i in range(size)]
            qubit_v_indices1 = [self._qregs.get_virtual_index(qreg1.name, i) for i in range(size)]

        for vi0, vi1 in zip(qubit_v_indices0, qubit_v_indices1):
            hi0 = self._qregs.convert_to_hardware_index(vi0)
            hi1 = self._qregs.convert_to_hardware_index(vi1)

            self._move_qubit_to_neighbor(hi0, hi1)

            new_hi0 = self._qregs.convert_to_hardware_index(vi0)
            new_hi1 = self._qregs.convert_to_hardware_index(vi1)

            self._engine.apply(CNOT(new_hi0, new_hi1))

    def _call_custom_unitary(self, args, env):
        """
        Calls a custom unitary gate (user-defined gate).

        Args:
            args (list[qiskit.qasm.node.Node]): The custom unitary gate info to be called.
            env (dict[str, qiskit.qasm.node.Node): The variable environment in which
                variable are evaluated to classical values or QASM register IDs.

        Raises:
            QASMError: Invalid QASM is used in the definition of the custom unitary gate.
        """

        gate_name = args[0].name
        gate_info = self._custom_unitaries[gate_name]

        inner_env = {}
        if len(gate_info) > 3:
            # If the gate receives some c-number arguments
            reg_arg_pos = 2

            for var, exp in zip(gate_info[1].children, args[1].children):
                inner_env[var.name] = qiskit.qasm.node.Real(exp.sym([env]))
        else:
            # If the gate receives no c-number arguments
            reg_arg_pos = 1

            # Gate receives at least one register argument.
            # (Because if not, what is the gate used for?)

        for reg_var, reg_id in zip(gate_info[reg_arg_pos].children, args[reg_arg_pos].children):
            if not env:
                inner_env[reg_var.name] = reg_id
            else:
                inner_env[reg_var.name] = env[reg_id.name]

        for stat in gate_info[reg_arg_pos+1].children:
            try:
                self._execute_statement(stat, inner_env)
            except QASMError as ex:
                raise QASMError(ex.get_op_code + ', used in ' + gate_name,
                                ex.get_reason())

    def _process_if_statement(self, args, env):
        """
        Evaluates "if statement".

        Args:
            args (list[qiskit.qasm.node.Node]): The custom unitary gate info to be called.
            env (dict[str, qiskit.qasm.node.Node): The variable environment in which
                variable are evaluated to classical values or QASM register IDs.

        Raises:
            QASMError: The classical registe is not found.
                (This would never be raised because the same check is performed at parsing phase).

        Note:
            The OpenQASM 2.0 allows "if statement" only at top level so the `env` may be unnecessary.
        """

        id = args[0].name
        value = args[1].value
        if id not in self._cregs:
            raise QASMError('if', 'Classical register {} not found'.format(id))

        if self._cregs.get_data(id) == value:
            self._execute_statement(args[2], env)

    def _move_qubit_to_neighbor(self, hardware_origin, hardware_target):
        """
        Moves a qubit to a neighboring site of another qubit.

        The move is achieved by sequential swapping of two neighboring qubits.

        Args:
            hardware_origin (int): The hardware qubit-index that one of its neighboring site
                will be replaced by `hardware_target`.
            hardware_target (int): The hardware qubit-index to be moved.
        """

        swap_list = self._engine.get_swap_path(hardware_origin, hardware_target)

        for i in range(len(swap_list)-1):
            self._engine.apply(Swap(swap_list[i], swap_list[i+1]))

        self._qregs.update_qubit_position(swap_list)
