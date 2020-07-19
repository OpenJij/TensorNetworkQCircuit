## This file isn't used currently.

import qiskit.qasm.node

def eval_expression(exp, env):
    if type(exp) in (qiskit.qasm.node.Int, qiskit.qasm.node.Real):
        return exp.sym()
    elif type(exp) == qiskit.qasm.node.Id:
        return env[exp.name]
    elif type(exp) == qiskit.qasm.node.Prefix:
        return eval_prefix(exp.children, env)
    elif type(exp) == qiskit.qasm.node.BinaryOp:
        return eval_binary_operator(exp.children, env)
    else:
        raise Exception("Unimplemented node")


def eval_prefix(exp, env):
    arg = eval_expression(exp[1], env)

    op = exp[0].operation() # operation() returns Python operator object
    return op(arg)


def eval_binary_operator(exp, env):
    arg1 = eval_expression(exp[1], env)
    arg2 = eval_expression(exp[2], env)

    op = exp[0].operation() # operation() returns Python operator object
    return op(arg1, arg2)
