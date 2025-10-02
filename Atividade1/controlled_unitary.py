from zyz_decomposition import zyz_decomposition


def controlled_unitary(unitaryMatrix):
    """
    Apply a controlled unitary operation to a 2x2 unitary matrix and return the circuit in OpenQASM.
    """

    # 1. Decompose the unitary matrix into ZYZ decomposition.
    angles = zyz_decomposition(unitaryMatrix)
    theta = angles["beta"]  # beta corresponds to theta in the decomposition
    phi = angles["alpha"]   # alpha corresponds to phi in the decomposition  
    lambdaa = angles["gamma"]  # gamma corresponds to lambda in the decomposition
    phase = angles["delta"]    # delta is the global phase
    A_op = {"Rz": lambdaa, "Ry": theta/2.0}
    B_op = {"Ry": -theta/2.0, "Rz": -(phi+lambdaa)/2.0}
    C_op = {"Rz": (phi - lambdaa)/2.0}
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

rz({A["Rz"]}) q[0];
ry({A["Ry"]}) q[0];
p({phase}) q[1];

cx q[1], q[0];

ry({B["Ry"]}) q[0];
rz({B["Rz"]}) q[0];

cx q[1], q[0];

rz({C["Rz"]}) q[0];


measure q -> c;
"""
    return circuit
