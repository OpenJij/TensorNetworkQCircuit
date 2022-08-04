class QASMError(Exception):
    def __init__(self, opcode, reason):
        super(QASMError, self).__init__(opcode + ': ' + reason)
        self._opcode = opcode
        self._reason = reason

    def get_op_code(self):
        return self.opcode

    def get_reason(self):
        return self.reason
