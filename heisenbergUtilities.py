import pandas as pd
import json
import numpy as np
from qiskit.quantum_info import SparsePauliOp, Operator, Pauli, Clifford
from qiskit.exceptions import QiskitError
from qiskit import QuantumCircuit

"""
This method loads a program from a QASM file. It was written and entirely belongs to the authors of SB-QOPS
"""
def load_program(name,path):
    qc = QuantumCircuit.from_qasm_file("{}/{}".format(path,name))
    qc.remove_final_measurements()
    if len(qc.clbits) > 0:
        for i in range(len(qc.clbits)):
            qc.measure(i, i)
    else:
        qc.measure_all()
    qc.remove_final_measurements()
    return qc.copy()

def evolve_pauli_exact(pauli_label, unitary):
    """Return exact Pauli expansion after conjugation"""
    op = SparsePauliOp(pauli_label).to_operator()
    evolved = unitary.adjoint().compose(op).compose(unitary)
    sp = SparsePauliOp.from_operator(evolved).simplify()

    return {
        p.to_label(): c for p, c in zip(sp.paulis, sp.coeffs)
    }

#Rules need verification
def give_transfer_rules():
    """
    Pauli transfer rules for Clifford gates.
    Format:
        {gate: {input_pauli: [(output_pauli, coefficient), ...]}}
    """

    transferRules = {

        # =========================
        # 1-qubit Clifford gates
        # =========================

        "h": {
            "I": [("I", 1)],
            "X": [("Z", 1)],
            "Y": [("Y", -1)],
            "Z": [("X", 1)],
        },

        "s": {  # phase gate (your "p")
            "I": [("I", 1)],
            "X": [("Y", 1)],
            "Y": [("X", -1)],
            "Z": [("Z", 1)],
        },

        "x": {
            "I": [("I", 1)],
            "X": [("X", 1)],
            "Y": [("Y", -1)],
            "Z": [("Z", -1)],
        },

        "y": {
            "I": [("I", 1)],
            "X": [("X", -1)],
            "Y": [("Y", 1)],
            "Z": [("Z", -1)],
        },

        "z": {
            "I": [("I", 1)],
            "X": [("X", -1)],
            "Y": [("Y", -1)],
            "Z": [("Z", 1)],
        },

        # =========================
        # 2-qubit Clifford gates
        # Qiskit ordering: control, target
        # =========================

        "cx": {
            "IX": [("IX", 1)],
            "IY": [("ZY", 1)],
            "IZ": [("ZZ", 1)],

            "XI": [("XX", 1)],
            "XX": [("XI", 1)],
            "XY": [("XZ", 1)],
            "XZ": [("XY", 1)],

            "YI": [("YX", 1)],
            "YX": [("YI", 1)],
            "YY": [("YZ", 1)],
            "YZ": [("YY", 1)],

            "ZI": [("ZI", 1)],
            "ZX": [("ZX", 1)],
            "ZY": [("ZY", 1)],
            "ZZ": [("ZZ", 1)],
        },

        "cz": {
            "II": [("II", 1)],
            "IZ": [("IZ", 1)],
            "ZI": [("ZI", 1)],
            "ZZ": [("ZZ", 1)],

            "XI": [("XI", 1)],
            "IX": [("IX", 1)],
            "XX": [("XX", -1)],

            "YI": [("YI", 1)],
            "IY": [("IY", 1)],
            "YY": [("YY", -1)],
        },
    }

    return transferRules

def try_transfer(gate_name, pauli_str, q_indices, num_qubits, transferRules):

    # extract only relevant qubits
    sub = extract_subpauli(pauli_str, q_indices, num_qubits)
    print(sub)

    gate_rules = transferRules[gate_name]

    if sub not in gate_rules:
        return {pauli_str: 1.0}

    outputs = {}

    for out_sub, coeff in gate_rules[sub]:

        full_out = insert_subpauli(
            list(pauli_str),
            q_indices,
            out_sub,
            num_qubits
        )
        print(full_out)

        outputs["".join(full_out)] = coeff

    return outputs

def is_instruction_clifford(inst):
    try:
        Clifford(inst)
        return True
    except QiskitError:
        return False
    
def extract_subpauli(pauli_str, q_indices, num_qubits):
    """Extract only gate-relevant qubits into a local string."""

    return "".join(
        pauli_str[num_qubits - 1 - q]
        for q in q_indices
    )


def insert_subpauli(full_pauli, q_indices, sub_pauli, num_qubits):

    full = list(full_pauli)

    for i, q in enumerate(q_indices):
        full[num_qubits - 1 - q] = sub_pauli[i]

    return "".join(full)

def create_single_test_dataframe(operation_list):
    node_list = []
    node = operation_list.head
    while node:
        node_list.append(node)
        node = node.next


    total_circuit_frame = pd.DataFrame({
        "total" : [node.count for node in node_list],
    })

    total_circuit_frame["total"] = total_circuit_frame["total"].astype(np.float64)

    return total_circuit_frame

def create_global_test_dataframe(operation_list):
    node_list = []
    node = operation_list.head
    while node:
        node_list.append(node)
        node = node.next

    global_circuit_frame = pd.DataFrame({}, columns = [node.value for node in node_list] + ["pass/fail"], dtype=np.float64)

    return global_circuit_frame

def format_test_case(test, first):
    #Cleaning up the test case string from the csv file to convert it into a dictionary
    test = test.replace("{", "")
    test = test.replace("]", "")
    test = test.replace("[", "")
    if not first:
        test = test[2:]
    test = "{" +test + "}"
    test = json.loads(test)

    return test

def add_counts_to_linked_list(operation_list, transition_graph, string_coeff):
    #Start at first gate in the inverse circuit
    checked_gate = operation_list.head.next
    idx = 1

    #For each gate in the circuit
    while checked_gate:

        #For each transition through that gate
        for edge in transition_graph[idx]["edges"]:
            
            #If the previous Pauli string evolved into a new Pauli string, 
            #add the probability of that new string occuring to that gate's count
            if edge["from"] != edge["to"]:
                checked_gate.count += edge["probability"]
        
        checked_gate.count *= abs(string_coeff)
        idx += 1
        checked_gate = checked_gate.next
    
    return operation_list

def append_to_analysis(testcase_analysis, operation_list, num_strings):
    node = operation_list.head
    gate_list = []

    while node:
        gate_list.append(node.count/num_strings)
        node = node.next
        
    gate_list.append("fail")
    testcase_analysis = pd.concat([testcase_analysis, pd.DataFrame([gate_list], columns=testcase_analysis.columns)], ignore_index=True)

    return testcase_analysis

def tarantula(testcase_analysis):
    num_fail_tests = len(testcase_analysis[testcase_analysis["pass/fail"] == "fail"])
    num_pass_tests = len(testcase_analysis[testcase_analysis["pass/fail"] == "pass"])
    fail_counts = testcase_analysis[testcase_analysis["pass/fail"] == "fail"].agg(["sum"]).drop(["pass/fail"], axis=1)
    pass_counts = testcase_analysis[testcase_analysis["pass/fail"] == "pass"].agg(["sum"]).drop(["pass/fail"], axis=1)

    tarantula_scores = (fail_counts/num_fail_tests)/((fail_counts/num_fail_tests)+(pass_counts/num_pass_tests))
    return tarantula_scores

def ochiai(testcase_analysis):
    pass_counts = testcase_analysis[testcase_analysis["pass/fail"] == "pass"].agg(["sum"]).drop(["pass/fail"], axis=1)
    fail_counts = testcase_analysis[testcase_analysis["pass/fail"] == "fail"].agg(["sum"]).drop(["pass/fail"], axis=1)
    num_fail_tests = len(testcase_analysis[testcase_analysis["pass/fail"] == "fail"])

    ochiai_scores = fail_counts/np.sqrt(num_fail_tests*(fail_counts+pass_counts))
    return ochiai_scores