from .exception import QASMError

class QuantumRegisters:
    def __init__(self, qubit_num):
        self._top = 0
        self._qubit_num = qubit_num
        self._qregs = {}
        self._index_map_virtual_to_hardware = [i for i in range(qubit_num)]


    def add(self, id, size):
        if id in self._qregs:
            raise QASMError('add_quantum_register', 'Multiple definition of quantum register {}'.format(id))

        if self._top + size > self._qubit_num:
            raise QASMError('add_quantum_register', 'Number of quantum registers exceeds total number of qubits')

        self._qregs[id] = QuantumRegisters.RegInfo(size, self._top)
        self._top += size


    def get_virtual_index(self, id, index):
        if index >= self._qregs[id].get_size():
            raise QASMError('get_hardware_index', 'Index exceeds register size')

        return self._qregs[id].get_start_index() + index


    def get_hardware_index(self, id, index):
        return self._index_map_virtual_to_hardware[self.get_virtual_index(id, index)]


    def convert_to_hardware_index(self, virtual_index):
        return self._index_map_virtual_to_hardware[virtual_index]


    def get_size(self, id):
        return self._qregs[id].get_size()


    def update_qubit_position(self, swap_list):
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
        if id in self._cregs:
            raise QASMError('add_classical_register', 'Multiple definition of classical register {}'.format(id))

        self._cregs[id] = ClassicalRegisters.RegInfo(size)

    def get(self, id, index):
        return self._cregs[id].get(index)


    def set(self, id, index, value):
        self._cregs[id].set(index, value)

    def get_data(self, id):
        return self._cregs[id].get_data()


    def get_size(self, id):
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
                index = self._size + index

            return (self._data >> index) & 1


        def set(self, index, value):
            if value == 0:
                self._data &= ~(1 << index)
            elif value == 1:
                self._data |= (1 << index)
            else:
                raise QASMError('set_classical_register_value', 'Classical registers cannot be set as {}'.format(value))


        def get_size(self):
            return self._size

