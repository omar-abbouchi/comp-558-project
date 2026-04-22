import cirq
import numpy as np

from circuits import circuit1, circuit2, circuit3

def total_variation_distance(dist1, dist2):
    return 0.5 * np.sum(np.abs(dist1 - dist2))

def normalize(value, max_value):
    return value / max_value if max_value > 0 else 0

def analyze_circuits(circuit1: cirq.Circuit, circuit2: cirq.Circuit, circuit3: cirq.Circuit):
    """
    Analyzes three given Cirq quantum circuits.

    Parameters:
        circuit1 (cirq.Circuit): The first quantum circuit to analyze.
        circuit2 (cirq.Circuit): The second quantum circuit to analyze.
        circuit3 (cirq.Circuit): The third quantum circuit to analyze.

    Returns:
        float: A single weighted score summarizing the performance of the three circuits.

    Raises:
        ValueError: If any circuit contains gates acting on more than two qubits, does not have exactly 8 qubits, or does not measure all 8 qubits.
    """
    circuits = [circuit1, circuit2, circuit3]
    target_distributions = [
        lambda n: np.power(np.arange(1, n+1), 0.2),  # Power-law distribution
        lambda n: np.exp(-0.5 * (np.arange(n) - n//2)**2 / (n//4)**2),  # Gaussian distribution
        lambda n: np.linspace(1, 0.1, n)  # Linearly decreasing distribution
    ]
    num_shots = 100000
    num_qubits = 8

    simulator = cirq.Simulator()
    scores = []

    max_tvd = 1  # TVD ranges from 0 to 1
    max_one_qubit = 100  # Assume a reasonable max for normalization
    max_two_qubit = 100  # Assume a reasonable max for normalization
    max_depth = 100  # Assume a reasonable max for normalization

    for i, circuit in enumerate(circuits):
        qubits = sorted(circuit.all_qubits())
        if len(qubits) != num_qubits:
            return 0

        # Check that all qubits are measured
        measured_qubits = set()
        for moment in circuit:
            for op in moment.operations:
                if isinstance(op.gate, cirq.MeasurementGate):
                    measured_qubits.update(op.qubits)
        if len(measured_qubits) != num_qubits:
            return 0

        one_qubit_gate_count = 0
        two_qubit_gate_count = 0

        for moment in circuit:
            for op in moment.operations:
                if isinstance(op.gate, cirq.MeasurementGate):
                    continue  # Skip measurement gates
                qubit_count = len(op.qubits)
                if qubit_count == 1:
                    one_qubit_gate_count += 1
                elif qubit_count == 2:
                    two_qubit_gate_count += 1
                else:
                    return 0

        circuit_depth = len(cirq.Circuit(circuit.all_operations()))

        # Simulate the circuit
        result = simulator.run(circuit, repetitions=num_shots)
        counts = result.histogram(key=list(result.measurements.keys())[0])

        # Convert counts to probabilities
        num_states = 2 ** num_qubits
        probabilities = np.zeros(num_states)
        for key, value in counts.items():
            probabilities[key] = value / num_shots

        # Compute TVD against the target distribution
        target_dist = target_distributions[i](num_states)
        target_dist /= np.sum(target_dist)  # Normalize
        tvd = total_variation_distance(probabilities, target_dist)

        # Normalize metrics
        norm_tvd = normalize(tvd, max_tvd)
        norm_one_qubit = normalize(one_qubit_gate_count, max_one_qubit)
        norm_two_qubit = normalize(two_qubit_gate_count, max_two_qubit)
        norm_depth = normalize(circuit_depth, max_depth)

        # Compute weighted score for this circuit
        circuit_score = 1 / (0.75 * norm_tvd + 0.05 * norm_one_qubit + 0.1 * norm_two_qubit + 0.1 * norm_depth)
        scores.append(circuit_score)

    # Average score across all circuits
    return sum(scores) / len(scores)

if __name__ == '__main__':

    score = analyze_circuits(circuit1, circuit2, circuit3)

    print('Score: ' + str(score))