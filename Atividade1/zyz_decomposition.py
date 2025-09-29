import numpy as np

def zyz_decomposition(unitaryMatrix):
    """
    Decompose a 2x2 unitary matrix into ZYZ decomposition.

    Input: 2x2 Unitary Matrix.
    Output: OpenQASM circuit decomposition using Ry, Rz and a global phase factor.  
    """

    # 1. Decomposing unitaryMatrix into theta, phi and psi.




    # 2. Implementing the decomposition into OpenQASM.
    circuit = f"""
OPENQASM 3.0;
include "stdgates.inc";

qubit[1] q;
bit[1] c;

// ZYZ Decomposition: U = Rz(phi) * Ry(theta) * Rz(psi)
rz({psi:.6f}) q[0];
ry({theta:.6f}) q[0];
rz({phi:.6f}) q[0];

measure q -> c;
"""
    
    return circuit
