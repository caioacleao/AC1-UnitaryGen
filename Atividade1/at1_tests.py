import numpy as np
from qiskit.quantum_info import Operator
from qiskit.synthesis import OneQubitEulerDecomposer
from qiskit.qasm3 import loads as qasm3_loads

from zyz_decomposition import zyz_decomposition, zyz_decomposition_circuit


"""
Testes para a decomposição ZYZ e síntese do operador controlado.

Os testes são realizados com 100 matrizes unitárias aleatórias.

A decomposição ZYZ e a síntese do operador controlado são comparados com a implementação do Qiskit.

Utilização:  

python -m pytest at1_tests.py -v -s

Caso o teste falhe, o parâmetro -s pode ser utilizado para mostrar os detalhes da comparação.

Requerimentos: pytest, qiskit, numpy e qiskit-qasm3.

pip install pytest, qiskit, qiskit_qasm3_import, numpy


"""

def unitary_from_angles(tetha, phi, lambdaa, phase):
    """Construct the 2x2 unitary from the ZYZ angles.
    
    Input: theta, phi, lambdaa, phase - ZYZ decomposition angles and global phase.
    Output: 2x2 Unitary Matrix, complex128 entries.
    
    The decomposition follows Qiskit convention: U = e^(i*phase) * Rz(phi) * Ry(theta) * Rz(lambda).
    """

    # Rz
    rz_phi = np.array(
        [[np.exp(-1j * phi / 2), 0], [0, np.exp(1j * phi / 2)]], dtype=np.complex128
    )
    
    # Ry
    ry_theta = np.array(
        [
            [np.cos(tetha / 2), -np.sin(tetha / 2)],
            [np.sin(tetha / 2), np.cos(tetha / 2)],
        ],
        dtype=np.complex128,
    )
    
    # Rz 
    rz_lambda = np.array(
        [[np.exp(-1j * lambdaa / 2), 0], [0, np.exp(1j * lambdaa / 2)]], dtype=np.complex128
    )

    # 
    unitaryMatrix = np.exp(1j * phase) * rz_phi @ ry_theta @ rz_lambda
    return unitaryMatrix


def random_unitary(rng):
    """Generate a random 2x2 unitary matrix using QR decomposition."""

    x = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    q, r = np.linalg.qr(x)
    d = np.diag(r)
    ph = d / np.abs(d)
    q = q * ph
    return q


def compare_to_qiskit(unitary, show_details=False):
    """Compare local ZYZ decomposition with Qiskit's implementation.
    
    Test that:
    1. Our decomposition correctly reconstructs the original matrix
    2. Qiskit's decomposition also correctly reconstructs the original matrix
    3. The OpenQASM circuit produces the correct unitary
    """

    local_angles = zyz_decomposition(unitary)
    theta, phi, lambdaa, phase = local_angles

    decomposer = OneQubitEulerDecomposer(basis="ZYZ")
    qiskit_theta, qiskit_phi, qiskit_lambda, qiskit_phase = decomposer.angles_and_phase(unitary)

    # Reconstruct from our angles
    local_u = unitary_from_angles(*local_angles)
    
    # Reconstruct from Qiskit angles
    # Qiskit returns (theta, phi, lambda, phase)
    qiskit_u = unitary_from_angles(qiskit_theta, qiskit_phi, qiskit_lambda, qiskit_phase)

    # Get unitary from OpenQASM circuit
    qasm = zyz_decomposition_circuit(*local_angles)
    qc = qasm3_loads(qasm)
    qc_no_measure = qc.remove_final_measurements(inplace=False)
    qasm_u = Operator(qc_no_measure).data

    # Check that our decomposition reconstructs the original
    our_reconstruction_ok = np.allclose(unitary, local_u, atol=1e-8)
    
    # Check that Qiskit's decomposition also reconstructs the original
    qiskit_reconstruction_ok = np.allclose(unitary, qiskit_u, atol=1e-8)
    
    # Check that the QASM circuit matches our reconstruction
    qasm_circuit_ok = np.allclose(local_u, qasm_u, atol=1e-8)

    # Show detailed comparison only if requested
    if show_details:
        print("\n" + "="*80)
        print("COMPARAÇÃO DETALHADA:")
        print("="*80)
        print("Matriz original:")
        print(unitary)
        print()
        
        print("Nossa decomposição:")
        print(f"  θ (theta)  = {theta}")
        print(f"  φ (phi)    = {phi}")
        print(f"  λ (lambda) = {lambdaa}")
        print(f"  γ (phase)  = {phase}")
        print()
        
        print("Decomposição Qiskit (theta, phi, lambda, phase):")
        print(f"  θ (theta)  = {qiskit_theta:.10f}")
        print(f"  φ (phi)    = {qiskit_phi:.10f}")
        print(f"  λ (lambda) = {qiskit_lambda:.10f}")
        print(f"  γ (phase)  = {qiskit_phase:.10f}")
        print()
        
        print("Nossa reconstrução:")
        print(local_u)
        print()
        
        print("Reconstrução Qiskit:")
        print(qiskit_u)
        print()
        
        print("Circuito QASM:")
        print(qasm_u)
        print()
        
        # Calculate differences
        diff_original_local = np.max(np.abs(unitary - local_u))
        diff_original_qiskit = np.max(np.abs(unitary - qiskit_u))
        diff_local_qiskit = np.max(np.abs(local_u - qiskit_u))
        diff_local_qasm = np.max(np.abs(local_u - qasm_u))
        
        print("Diferenças:")
        print(f"  Original vs Nossa:      {diff_original_local:.2e}")
        print(f"  Original vs Qiskit:     {diff_original_qiskit:.2e}")
        print(f"  Nossa vs Qiskit:        {diff_local_qiskit:.2e}")
        print(f"  Nossa vs QASM:          {diff_local_qasm:.2e}")
        print()
        
        print("Resultados:")
        print(f"  Nossa reconstrução OK:     {our_reconstruction_ok}")
        print(f"  Qiskit reconstrução OK:     {qiskit_reconstruction_ok}")
        print(f"  QASM circuito OK:           {qasm_circuit_ok}")
        print("="*80)

    return qiskit_reconstruction_ok and our_reconstruction_ok and qasm_circuit_ok

def test_zyz_decomposition_matches_qiskit():
    """Test 100 random unitaries against Qiskit's decomposer."""

    rng = np.random.default_rng()
    failed_matrices = []
    
    for i in range(100):
        u = random_unitary(rng)
        result = compare_to_qiskit(u, show_details=False)
        
        if not result:
            failed_matrices.append((i, u))
    
    if failed_matrices:
        print(f"\n Teste falhou em {len(failed_matrices)} matrizes. Mostrando detalhes da primeira falha (matriz {failed_matrices[0][0]+1}):")
        compare_to_qiskit(failed_matrices[0][1], show_details=True)
        assert False, f"Teste falhou em {len(failed_matrices)} matrizes. Primeira falha na matriz {failed_matrices[0][0]+1}."
    else:
        print(f"{100 - len(failed_matrices)}/100 testes passaram")
