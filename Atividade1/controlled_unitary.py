from zyz_decomposition import zyz_decomposition


def controlled_unitary(unitaryMatrix):
    """
    Apply a controlled unitary operation to a 2x2 unitary matrix and return the circuit in OpenQASM.
    """

    # 1. Decompose the unitary matrix into ZYZ decomposition.
    theta, phi, lambdaa, phase = zyz_decomposition(unitaryMatrix)

    A_op = {"Rz": lambdaa, "Ry": theta / 2.0}
    B_op = {"Ry": -theta / 2.0, "Rz": -(phi + lambdaa) / 2.0}
    C_op = {"Rz": (phi - lambdaa) / 2.0}
    return A_op, B_op, C_op, phase


def controlled_unitary_circuit(unitaryMatrix):
    """
    Implement the controlled unitary operation in OpenQASM.
    """

    A, B, C, phase = controlled_unitary(unitaryMatrix)

    circuit = f"""
OPENQASM 3.0;
include "stdgates.inc";

qubit[2] q;
bit[2] c;

rz({C["Rz"]:.6f}) q[1];

cx q[0], q[1];

rz({B["Rz"]:.6f}) q[1];
ry({B["Ry"]:.6f}) q[1];

cx q[0], q[1];

ry({A["Ry"]:.6f}) q[1];
rz({A["Rz"]:.6f}) q[1];
p({phase:.6f}) q[0];

measure q -> c;
"""
    return circuit
