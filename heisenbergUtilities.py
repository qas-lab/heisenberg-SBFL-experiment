import pandas as pd
import json
import numpy as np
from qiskit.quantum_info import SparsePauliOp, Operator, Pauli, Clifford
from qiskit.exceptions import QiskitError
from qiskit import QuantumCircuit
from pauli_prop import propagate_through_operator

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

def evolve_pauli_exact(pauli_label, gate):
    """Return exact Pauli expansion after conjugation"""
    evolution = propagate_through_operator(pauli_label, gate, atol=1e-4, frame='h', max_terms = 80, search_step=3)

    return {
        p.to_label(): c for p, c in zip(evolution.paulis, evolution.coeffs)
    }

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

def append_to_analysis(testcase_analysis, operation_list, num_strings, pass_fail):
    node = operation_list.head
    gate_list = []

    while node:
        gate_list.append(node.count/num_strings)
        node = node.next
    
    if pass_fail == "fail":
        gate_list.append("fail")
    else:
        gate_list.append("pass")
        
    testcase_analysis = pd.concat([testcase_analysis, pd.DataFrame([gate_list], columns=testcase_analysis.columns)], ignore_index=True)

    return testcase_analysis

def tarantula(testcase_analysis):
    num_fail_tests = len(testcase_analysis[testcase_analysis["pass/fail"] == "fail"])
    num_pass_tests = len(testcase_analysis[testcase_analysis["pass/fail"] == "pass"])
    fail_counts = testcase_analysis[testcase_analysis["pass/fail"] == "fail"].agg(["sum"]).drop(["pass/fail"], axis=1)
    pass_counts = testcase_analysis[testcase_analysis["pass/fail"] == "pass"].agg(["sum"]).drop(["pass/fail"], axis=1)

    tarantula_scores = (fail_counts/num_fail_tests)/((fail_counts/num_fail_tests)+(pass_counts/num_pass_tests))
    tarantula_scores = tarantula_scores[tarantula_scores.iloc[0].sort_values(ascending=False).index]
    return tarantula_scores

def ochiai(testcase_analysis):
    pass_counts = testcase_analysis[testcase_analysis["pass/fail"] == "pass"].agg(["sum"]).drop(["pass/fail"], axis=1)
    fail_counts = testcase_analysis[testcase_analysis["pass/fail"] == "fail"].agg(["sum"]).drop(["pass/fail"], axis=1)
    num_fail_tests = len(testcase_analysis[testcase_analysis["pass/fail"] == "fail"])

    ochiai_scores = fail_counts/np.sqrt(num_fail_tests*(fail_counts+pass_counts))
    ochiai_scores = ochiai_scores[ochiai_scores.iloc[0].sort_values(ascending=False).index]
    return ochiai_scores