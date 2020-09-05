from .exception import QASMError

class QuantumRegisters:
    """
    The class to store infomation of defined quantum registers.

    This class stores quantum registers as a "virtual" one dimensional array,
    e.g. if one declares "qreg q[3];" and then "qreg r[4];", the quantum register "q"
    uses the virtual index from 0 to 2 and the "r" uses from 3 to 6.
    As well as this example, the virtual indices are consumed in ascending order.
    A virtual index corresponds to a hardware index at a certain moment,
    but the correspondence changes dynamically.

    Attributes:
        _top (int): The minimum index not being used, e.g. in the above example
            this variable is 7.
        _qubit_num(int): The number of hardware qubits.
        _qregs (dict[str, QuantumRegisters.RegInfo]): The dictionary to convert a quantum register name
            to its info.
        _index_map_virtual_to_hardware (list[int]): The list which maps a virtual index to a hardware index.
    """

    def __init__(self, qubit_num):
        """
        Constructs a object to store quantum register information.

        Args:
            qubit_num (int): The number of hardware qubits,
                which limits the number of total qubits used by registers.
        """

        self._top = 0
        self._qubit_num = qubit_num
        self._qregs = {}
        self._index_map_virtual_to_hardware = [i for i in range(qubit_num)]

    def add(self, id, size):
        """
        Add a quantum register.

        Args:
            id (str): The name of the register.
            size (int): The size of the register.

        Raises:
            QASMError: The name is already used.
                Or the number of qubits being registered exceeds the number of hardware qubits.
        """

        if id in self._qregs:
            raise QASMError('add_quantum_register',
                            'Multiple definition of quantum register {}'.format(id))

        if self._top + size > self._qubit_num:
            raise QASMError('add_quantum_register',
                            'Number of quantum registers exceeds total number of qubits')

        self._qregs[id] = QuantumRegisters.RegInfo(size, self._top)
        self._top += size

    def get_virtual_index(self, id, index):
        """
        Returns the virtual index corresponding to the given quantum register.

        More about the virtual index, see the class documentation.

        Args:
            id (str): The name of the register.
            index (int): The index of the register.

        Returns:
            int: The virtual index.

        Raises:
            QASMError: The index exceeds the size of the register.
        """

        if index >= self._qregs[id].get_size():
            raise QASMError('get_hardware_index',
                            'Index exceeds register size')

        return self._qregs[id].get_start_index() + index

    def get_hardware_index(self, id, index):
        """
        Returns the hardware index corresponding to the given quantum register.

        More about the hardware index, see the class documentation.

        Args:
            id (str): The name of the register.
            index (int): The index of the register.

        Returns:
            int: The hardware index.

        Raises:
            QASMError: The index exceeds the size of the register.
        """

        return self._index_map_virtual_to_hardware[self.get_virtual_index(id, index)]

    def convert_to_hardware_index(self, virtual_index):
        """
        Returns the hardware index corresponding to the given virtual index.

        More about the virtual/hardware index, see the class documentation.

        Args:
            virtual_index (int): The virtual index.

        Returns:
            int: The hardware index.
        """

        return self._index_map_virtual_to_hardware[virtual_index]

    def get_size(self, id):
        """
        Returns the size of the quantum register.

        Args:
            id (str): The name of the register.

        Returns:
            int: The register size.
        """

        return self._qregs[id].get_size()

    def update_qubit_position(self, swap_list):
        """
        Updates qubit position by swapping qubits.

        This method change the correspondence between the virtual indices and
        the hardware indices, swapping ith and (i+1)th hardware index from left to right.
        As an example, if that the initial correspondence map
        (`_index_map_virtual_to_hardware`) is [0 1 2 3 4] and [4, 3, 1, 0] is given as an argument,
        The process is as follows

        1. Swaps 4 and 3, thus the correspondence map becomes [0 1 2 4 3]
        2. Swaps 3 and 1, thus the correspondence map becomes [0 3 2 4 1]
        3. Swaps 1 and 0, thus the correspondence map becomes [1 3 2 4 0]

        Args:
            swap_list (list[int]): The list of hardware indices.

        Returns:
            int: The register size.
        """

        for i in range(len(swap_list)-1):
            virtual_index0 = self._index_map_virtual_to_hardware.index(swap_list[i])
            virtual_index1 = self._index_map_virtual_to_hardware.index(swap_list[i+1])

            # swap items
            tmp = self._index_map_virtual_to_hardware[virtual_index0]
            self._index_map_virtual_to_hardware[virtual_index0] = self._index_map_virtual_to_hardware[virtual_index1]
            self._index_map_virtual_to_hardware[virtual_index1] = tmp

    def __contains__(self, item):
        return item in self._qregs

    class RegInfo:
        def __init__(self, size, start_index):
            self._size = size
            self._start_index = start_index

        def get_size(self):
            return self._size

        def get_start_index(self):
            return self._start_index

        def __repr__(self):
            return "RegInfo(start={}, size={})".format(self.start, self.size)



class ClassicalRegisters:
    def __init__(self):
        self._cregs = {}

    def add(self, id, size):
        """
        Add a classical register.

        Args:
            id (str): The name of the register.
            size (int): The size of the register.

        Raises:
            QASMError: The name is already used.
        """

        if id in self._cregs:
            raise QASMError('add_classical_register',
                            'Multiple definition of classical register {}'.format(id))

        self._cregs[id] = ClassicalRegisters.RegInfo(size)

    def get(self, id, index):
        """
        Returns the classical bit value, 0 or 1.

        Args:
            id (str): The name of the register.
            index (int): The index of the register.

        Returns:
            int: The classical bit value, 0 or 1.

        Raise:
            QASMError: The index exceeds the number of bits.
        """

        return self._cregs[id].get(index)

    def set(self, id, index, value):
        """
        Returns the classical bit value.

        Args:
            id (str): The name of the register.
            index (int): The index of the register.
            value (int): The value to be set, which should be 0 or 1.

        Raise:
            QASMError: The index exceeds the number of bits.
                Or The given value is neither 0 nor 1.
        """

        self._cregs[id].set(index, value)

    def get_data(self, id):
        """
        Returns the classical register value.

        Args:
            id (str): The name of the register.

        Returns:
            int: The classical register value,
                 which is in the range [0, 2**N-1] (N is the number of classical bits).
        """

        return self._cregs[id].get_data()

    def get_size(self, id):
        """
        Returns the size of the classical register.

        Args:
            id (str): The name of the register.

        Returns:
            int: The size of the classical register.
        """

        return self._cregs[id].get_size()

    def __contains__(self, item):
        return item in self._cregs

    class RegInfo:
        def __init__(self, size):
            self._size = size
            self._data = 0

        def get_data(self):
            return self._data

        def get(self, index):
            if index < 0:
                index = self._size + index  # negative index access (backward).

            return (self._data >> index) & 1

        def set(self, index, value):
            if value == 0:
                self._data &= ~(1 << index)
            elif value == 1:
                self._data |= (1 << index)
            else:
                raise QASMError('set_classical_register_value',
                                'Classical registers cannot be set as {}'.format(value))

        def get_size(self):
            return self._size
