import numpy as np

def zyz_decomposition(unitaryMatrix):
    
    """
    Decompose a 2x2 unitary matrix into ZYZ decomposition.

    Input: 2x2 Unitary Matrix, complex128 entries.
    Output: Rz, Ry, Rz and global phase factor.
    
    The decomposition follows Qiskit convention: U = e^(i*phase) * Rz(phi) * Ry(theta) * Rz(lambda)
    """

    # Extract global phase (determinant phase)
    det = np.linalg.det(unitaryMatrix)
    phase = np.angle(det) / 2
    
    # Remove global phase to get SU(2) matrix
    U_prime = unitaryMatrix * np.exp(-1j * phase)
    
    # Extract matrix elements
    u00 = U_prime[0, 0]  # The other 2 elements are not used in the decomposition.
    u10 = U_prime[1, 0]
    
    eps = 1e-10
    
    # Compute beta using |u00| = cos(beta/2)
    theta = 2 * np.arccos(np.abs(u00))
    
    phase_sum = -2 * np.angle(u00 / np.cos(theta / 2))
    phase_diff = 2 * np.angle(u10 / np.sin(theta / 2))
    
    phi = (phase_sum + phase_diff) / 2
    lambdaa = (phase_sum - phase_diff) / 2

    return (theta, phi, lambdaa, phase)


def zyz_decomposition_circuit(theta, phi, lambdaa, phase):
    """
    Implement the ZYZ decomposition into OpenQASM.

    Input: theta, phi, lambdaa, phase - ZYZ decomposition angles and global phase.
    Output: OpenQASM circuit decomposition using Ry, Rz and a global phase factor.  
    """

    circuit = f"""
OPENQASM 3.0;
include "stdgates.inc";

qubit[1] q;
bit[1] c;

// ZYZ Decomposition: U = e^(i*phase) * Rz(phi) * Ry(theta) * Rz(lambda).
// In OpenQASM, gates are applied from right to left


gphase({phase:.6f});
rz({lambdaa:.6f}) q[0];
ry({theta:.6f}) q[0];
rz({phi:.6f}) q[0];

measure q -> c;
"""
    return circuit
