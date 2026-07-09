import pandas as pd
import json
import numpy as np
from qiskit.quantum_info import SparsePauliOp, Operator, Pauli, Clifford
from qiskit.exceptions import QiskitError
from qiskit import QuantumCircuit, qasm3, qasm2
from pauli_prop import propagate_through_operator

"""
This method loads a program from a QASM file. It has only been modified to parse QASM3 files and otherwise belongs to the authors of SB-QOPS.

INPUTS:
    name (string): The name of the program file to be loaded. Must be QASM file.

    path (string): The filepath to get to the program from working directory.

OUTPUTS:
    qc (QuantumCircuit): A copy of the quantum circuit with measurements removed.
"""
def load_program(name,path):
    try:
        qc = qasm3.load("{}/{}".format(path,name))
        qc.remove_final_measurements()
        if len(qc.clbits) > 0:
            for i in range(len(qc.clbits)):
                qc.measure(i, i)
        else:
            qc.measure_all()
        qc.remove_final_measurements()
        return qc.copy()
    except:
        try:
            qc = qasm2.load("{}/{}".format(path,name), custom_instructions=qasm2.LEGACY_CUSTOM_INSTRUCTIONS)
            qc.remove_final_measurements()
            if len(qc.clbits) > 0:
                for i in range(len(qc.clbits)):
                    qc.measure(i, i)
            else:
                qc.measure_all()
            qc.remove_final_measurements()
            return qc.copy()
        except:
            raise Exception("File open error")

def construct_list(list, circuit, inverse):
    if inverse:
        depth = 0

        #Append each gate to the linked list
        list.append("Initial", depth)
        for instruction in circuit.data:
            depth += 1
            list.append(instruction.name + " " + str((circuit.size()-depth)), depth)
    else:
        depth = 0
        for instruction in circuit.data:
            list.append(instruction.name + " " + str(depth), depth)
            depth += 1

    return list

"""
This method utilizes the pauli-prop package to perform exact evolution over the desired gate.

INPUTS:
    pauli_label (SparsePauliOp): An operator representing the Pauli string that will evolve over the gate.

    gate (SparsePauliOp): The quantum gate over which the Pauli string will propagate. 

OUTPUTS:
    A dictionary object with the following form:
    {
    pauli1 : coeff1,
    pauli2 : coeff2,
    ...
    }
"""
def evolve_pauli_exact(pauli_label, gate, tolerance = 1e-4, terms = None, search_step = None):
    """Return exact Pauli expansion after conjugation"""
    evolution = propagate_through_operator(pauli_label, gate, atol=tolerance, frame='h', max_terms=terms, search_step=search_step)

    return {
        p.to_label(): c for p, c in zip(evolution.paulis, evolution.coeffs)
    }

"""
This method creates an overall pandas dataframe with columns that represent each gate in the quantum circuit as well as an indicator
of if the associated test/pauli was a passing or failing test

INPUTS:
    operation_list (LinkedList): A linked list structure of all the gates in the quantum circuit, including gate names, depth, and 
        probability counts.

OUTPUTS:
    global_circuit_frame (DataFrame): A pandas dataframe whose columns are gate labels + depth and a column for pass or fail indication
"""
def create_global_test_dataframe(operation_list):
    node_list = []
    node = operation_list.head
    while node:
        node_list.append(node)
        node = node.next

    global_circuit_frame = pd.DataFrame({}, columns = [node.value for node in node_list] + ["pass/fail"], dtype=np.float64)

    return global_circuit_frame

"""
This method formats the raw test cases into a normalized format that we can read from. 

INPUTS:
    test (string): An unfiltered string holding Pauli strings and coefficients for a test case.

    first (bool): A boolean statement specifying if this is the first test in the testcase batch. Used to offset characters that only 
        appear after the first test case.

OUTPUTS:
    test (string): The properly formatted version of the test input string
"""
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

"""
This method looks for potentially empty test cases and removes them so that it doesn't break our work flow or add faulty tests to the overall
results.

INPUTS:
    tests (List): A python list of unformatted test case strings.

OUTPUTS:
    tests (List): The original input list, but with any null entries removed.
"""
def remove_null_tests(tests):
    for raw_idx, raw in enumerate(tests):
        if raw == ']':
            tests.pop(raw_idx)

    return tests

"""

"""
def pauli_similarity(pauli_a, pauli_b):
    matches = sum(a == b for a, b in zip(pauli_a, pauli_b))
    return matches / len(pauli_a)

"""
This method adds the probability counts of any evolution that changed the Pauli string to the gate in the Linked List. 

INPUTS:
    operation_list (LinkedList): A Linked List of gates from the circuit that includes gate name, circuit depth, and a running count of the
        probability that it will change the Pauli string at a given moment

    transition_graph (Dict): A dictionary containing the change information regarding an evolution step from Pauli propagation

    string_coeff (Float): The coefficient from the test case indicating how to weight the final Pauli. We use this to experimentally
        weight the probability counts by however the test case weighs the Pauli we're analyzing

OUTPUTS:
    operation_list (LinkedList): The input Linked List with the counts updated for each gate 
"""
def add_counts_to_linked_list(operation_list, transition_graph, string_coeff, lambda_phase, lambda_change, target_pauli):
    #Start at first gate in the inverse circuit
    checked_gate = operation_list.head.next
    idx = 1

    #For each gate in the circuit
    while checked_gate:

        #For each transition through that gate
        for edge in transition_graph[idx]["edges"]:
            #-----------------------------------------------------------------------
            # score = 0

            # #If the previous Pauli string evolved into a new Pauli string, 
            # #add the probability of that new string occuring to that gate's count
            # if edge["from"] != edge["to"]:
            #     score += 1

            # score += edge["phase"] / np.pi

            # checked_gate.count += score * edge["probability"]
            #-----------------------------------------------------------------------

            similarity_difference = edge["to_similarity"] - edge["from_similarity"]

            distance = lambda_change * int(edge["from"]!=edge["to"]) + lambda_phase * abs(edge["phase"]) / np.pi

            checked_gate.count += (distance + similarity_difference) * edge["probability"]

        idx += 1
        checked_gate = checked_gate.next
    
    return operation_list

"""
This method appends a resulting Linked List with counts for all the gates to an overall analysis dataframe to use in SBFL

INPUTS:
    testcase_analysis (DataFrame): A pandas dataframe that is recording all of our test cases including gate counts and pass or fail

    operation_list (LinkedList): A Linked List containing counts for all gates after a test of Pauli evolutions

    num_strings (Int): The number of Pauli strings the specific test case utilized. We use this to weight the counts further,
        as to not bias our results towards test cases that required more Pauli strings to specify
    
    pass_fail (String): A string indicating if the test case we're appending was a passing or a failing test

OUTPUTS:
    testcase_analysis (DataFrame): An updated pandas DataFrame with the results of a test from the test case batch appended to it
"""
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

"""
This method is the implementation of the SBFL Tarantula algorithm, fitted to work with our data format.

INPUTS:
    testcase_analysis (DataFrame): A pandas DataFrame with the counts for all gates across all tests from our test cases.

OUTPUTS:
    tarantula_scores (DataFrame): A pandas DataFrame with columns ordered from highest to lowest suspiciousness scores based on the 
        Tarantula algorithm.
"""
def tarantula(testcase_analysis):
    num_fail_tests = len(testcase_analysis[testcase_analysis["pass/fail"] == "fail"])
    num_pass_tests = len(testcase_analysis[testcase_analysis["pass/fail"] == "pass"])
    fail_counts = testcase_analysis[testcase_analysis["pass/fail"] == "fail"].agg(["sum"]).drop(["pass/fail"], axis=1)
    pass_counts = testcase_analysis[testcase_analysis["pass/fail"] == "pass"].agg(["sum"]).drop(["pass/fail"], axis=1)

    tarantula_scores = (fail_counts/num_fail_tests)/((fail_counts/num_fail_tests)+(pass_counts/num_pass_tests))
    tarantula_scores = tarantula_scores[tarantula_scores.iloc[0].sort_values(ascending=False).index]
    return tarantula_scores

# """
# This method is the implementation of the SBFL Ochiai algorithm, fitted to work with our data format.

# INPUTS:
#     testcase_analysis (DataFrame): A pandas DataFrame with the counts for all gates across all tests from our test cases.

# OUTPUTS:
#     ochiai_scores (DataFrame): A pandas DataFrame with columns ordered from highest to lowest suspiciousness scores based on the Ochiai
#         algorithm.

# NOTE: I have discovered that Ochiai is not feasible to implement. It requires a differentiation between failing tests that involve a gate and failing
# tests that do not involve a gate. (Also this current implementation is incorrect)
# """
# def ochiai(testcase_analysis):
#     pass_counts = testcase_analysis[testcase_analysis["pass/fail"] == "pass"].agg(["sum"]).drop(["pass/fail"], axis=1)
#     fail_counts = testcase_analysis[testcase_analysis["pass/fail"] == "fail"].agg(["sum"]).drop(["pass/fail"], axis=1)
#     num_fail_tests = len(testcase_analysis[testcase_analysis["pass/fail"] == "fail"])

#     ochiai_scores = fail_counts/np.sqrt(num_fail_tests*(fail_counts+pass_counts))
#     ochiai_scores = ochiai_scores[ochiai_scores.iloc[0].sort_values(ascending=False).index]
#     return ochiai_scores

"""
This method is the implementation of the SBFL Barinel algorithm, fitted to work with our data format.

INPUTS:
    testcase_analysis (DataFrame): A pandas DataFrame with the counts for all gates across all tests from our test cases.

OUTPUTS:
    barinel_scores (DataFrame): A pandas DataFrame with columns ordered from highest to lowest suspiciousness scores based on the Barinel
        algorithm.

"""
def barinel(testcase_analysis):
    pass_counts = testcase_analysis[testcase_analysis["pass/fail"] == "pass"].agg(["sum"]).drop(["pass/fail"], axis=1)
    fail_counts = testcase_analysis[testcase_analysis["pass/fail"] == "fail"].agg(["sum"]).drop(["pass/fail"], axis=1)

    barinel_scores = 1 - ((pass_counts)/(fail_counts + pass_counts))
    barinel_scores = barinel_scores[barinel_scores.iloc[0].sort_values(ascending=False).index]
    return barinel_scores

def custom_sbfl(testcase_analysis):
    pass_counts = testcase_analysis[testcase_analysis["pass/fail"] == "pass"].agg(["sum"]).drop(["pass/fail"], axis=1)
    fail_counts = testcase_analysis[testcase_analysis["pass/fail"] == "fail"].agg(["sum"]).drop(["pass/fail"], axis=1)

    #TODO: Implement a custom SBFL algorithm here

"""
This method compares the original quantum circuit to the mutated one, and locates the depth where an added gate occurs.

INPUTS:
    forward_list (LinkedList): A Linked List containing the data regarding our NOT INVERSE mutant circuit

    correct_list (LinkedList): A linked list containing the data regarding the NOT INVERSE correct circuit

OUTPUTS:
    depth (Int): The depth of the erroneous gate in the mutant circuit.

NOTE: This method only works for mutants that add a gate to the circuit. Functionality for mutants that remove or
    replace a gate is for future work, as SBFL works best for mutants with extra gates
"""
def find_erroneous_gate(forward_list, correct_list):
    correct_head = correct_list.head
    mutant_head = forward_list.head
    while mutant_head:
        if correct_head is None:
            return mutant_head.depth
        if correct_head.value == mutant_head.value:
            correct_head = correct_head.next
            mutant_head = mutant_head.next
        else:
            return mutant_head.depth
    
    return 0