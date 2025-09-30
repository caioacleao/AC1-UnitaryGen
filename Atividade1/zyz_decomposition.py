import numpy as np

def zyz_decomposition(unitaryMatrix):
    
    """
    Decompose a 2x2 unitary matrix into ZYZ decomposition.

    Input: 2x2 Unitary Matrix, complex128 entries.
    Output: Rz, Ry, Rz and global phase factor.
    
    The decomposition follows Qiskit convention: U = e^(i*delta) * Rz(alpha) * Ry(beta) * Rz(gamma)
    """

    # Extract global phase (determinant phase)
    det = np.linalg.det(unitaryMatrix)
    delta = np.angle(det) / 2
    
    # Remove global phase to get SU(2) matrix
    U_prime = unitaryMatrix * np.exp(-1j * delta)
    
    # Extract matrix elements
    u00 = U_prime[0, 0]  # The other 2 elements are not used in the decomposition.
    u10 = U_prime[1, 0]
    
    eps = 1e-10
    
    # Compute beta using |u00| = cos(beta/2)
    beta = 2 * np.arccos(np.clip(np.abs(u00), 0, 1))
    
    # Handle special cases
    if np.abs(beta) < eps:
        # beta ≈ 0, matrix is approximately Rz(alpha+gamma)
        # U ≈ Rz(alpha+gamma), set gamma = 0
        alpha = 2 * np.angle(u00)
        gamma = 0.0
    elif np.abs(beta - np.pi) < eps:
        # beta ≈ π, matrix is approximately Rz(alpha-gamma) * Ry(π)
        # Set gamma = 0
        alpha = 2 * np.angle(u10)
        gamma = 0.0
    else:
        # u00 = cos(beta/2) * e^(-i(alpha+gamma)/2)
        # u01 = -sin(beta/2) * e^(-i(alpha-gamma)/2)
        # u10 = sin(beta/2) * e^(i(alpha-gamma)/2)
        # u11 = cos(beta/2) * e^(i(alpha+gamma)/2)
        
        # Extract phases
        # alpha + gamma = -2 * angle(u00 / cos(beta/2))
        # alpha - gamma = 2 * angle(u10 / sin(beta/2))
        
        phase_sum = -2 * np.angle(u00 / np.cos(beta / 2))
        phase_diff = 2 * np.angle(u10 / np.sin(beta / 2))
        
        alpha = (phase_sum + phase_diff) / 2
        gamma = (phase_sum - phase_diff) / 2

    return dict(alpha=alpha, beta=beta, gamma=gamma, delta=delta)


def zyz_decomposition_circuit(alpha, beta, gamma, delta):
    """
    Implement the ZYZ decomposition into OpenQASM.

    Input: alpha, beta, gamma, delta - ZYZ decomposition angles and global phase.
    Output: OpenQASM circuit decomposition using Ry, Rz and a global phase factor.  
    """

    circuit = f"""
OPENQASM 3.0;
include "stdgates.inc";

qubit[1] q;
bit[1] c;

// ZYZ Decomposition: U = e^(i*delta) * Rz(alpha) * Ry(beta) * Rz(gamma).
// In OpenQASM, gates are applied from right to left


gphase({delta:.6f});
rz({gamma:.6f}) q[0];
ry({beta:.6f}) q[0];
rz({alpha:.6f}) q[0];

measure q -> c;
"""
    return circuit