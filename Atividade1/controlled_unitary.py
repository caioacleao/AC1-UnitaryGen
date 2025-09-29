import numpy as np
from zyz_decomposition import zyz_decomposition


def controlled_unitary(unitaryMatrix):
    """
    Apply a controlled unitary operation to a 2x2 unitary matrix and return the circuit in OpenQASM.
    """

    # 1. Decompose the unitary matrix into ZYZ decomposition.
    theta, phi, psi = zyz_decomposition(unitaryMatrix)

    
    # 2. Implement the decomposition into OpenQASM.
    circuit = f"""
OPENQASM 3.0;
include "stdgates.inc";

qubit[2] q;
bit[2] c;

// Implementing the controlled unitary operation.
    return circuit
