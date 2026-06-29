import random
import threading
from datetime import datetime

import numpy as np
import pandas as pd
from mthree.utils import counts_to_vector, vector_to_quasiprobs
from pytket._tket.partition import MeasurementSetup, measurement_reduction, PauliPartitionStrat, MeasurementBitMap
from pytket.backends.backendresult import BackendResult

from qiskit import generate_preset_pass_manager, transpile
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp
from qiskit_aer import AerSimulator
from qiskit_ibm_runtime import Session, EstimatorOptions
from qiskit_ibm_runtime.estimator import EstimatorV2 as Estimator

from MyCustomES import customES, SearchSpace
from QOPS_test import get_random_circuit, get_Z_family_values, get_random_Z_family, get_compact_program_specification_Z
from mealpy import FloatVar, GA, HC, SA, TS, IntegerVar, ES, EP
from mealpy import Problem as MealProblem
import time
from qiskit_ibm_runtime import QiskitRuntimeService, EstimatorV2 as IBMestimator

class Circuit_Tester:
    def __init__(self,CUT:QuantumCircuit,simulator_type="statevector",shots=20000):
        """
        :param CUT: A qiskit circuit without measurements
        """
        self.compact_program_specification = None
        self.CUT = CUT
        self.aer_sim = AerSimulator(method=simulator_type, device='GPU', blocking_enable=True, batched_shots_gpu=True)
        self.pass_manager = generate_preset_pass_manager(backend=self.aer_sim, optimization_level=3)
        self.shots = shots

    def set_applicable_families_Z(self,compact_program_specification:dict):
        """
        :param compact_program_specification: dict of paulistring with their measurement results
        """
        self.compact_program_specification = compact_program_specification
        self.applicable_families = []
        zfam = "Z"*(2**self.CUT.num_qubits)
        for key,value in compact_program_specification.items():
            if key in zfam:
                self.applicable_families.append((0, value))

    def random_test_case_Z(self):
        pauli_dict = {}
        fam = random.choice(self.applicable_families)
        selected_pauli = get_random_Z_family(self.CUT.num_qubits)
        for s in selected_pauli:
            pauli_dict[s] = np.random.random()
        return {"test_case":pauli_dict,"family_index":fam[0],"M":fam[-1]}


    def get_test_case_theoretics_Z(self,test_case:dict):
        return {"test_case": test_case["test_case"], "M": test_case["M"]}

    def get_theoretical_exp_from_testcase(self, test_case:dict):

        def pauli_expectation_from_counts(counts, pauli_string):
            pauli_string = pauli_string.upper()
            total_counts = sum(counts.values())
            exp_val = 0.0
            for bitstring, count in counts.items():
                v = 1
                for i, p in enumerate(pauli_string):
                    if p == 'Z':
                        v *= 1 if bitstring[i] == '0' else -1
                exp_val += v * count
            return exp_val / total_counts

        # Compute expectation value of the Hamiltonian
        counts = test_case["M"]
        expectation = sum(
            coeff * pauli_expectation_from_counts(counts, pauli_str)
            for pauli_str, coeff in test_case['test_case'].items()
        )
        return expectation

    def get_theoretical_exp_from_test_case_M3(self, test_case:dict):

        counts = test_case["M"]
        quasi = vector_to_quasiprobs(counts_to_vector(counts),counts)
        expectation = sum([v*quasi.expval(exp_ops=x) for x,v in test_case["test_case"].items()])
        return expectation

    def execute_test_cases(self,test_cases:list):
        pubs = []
        for test_case in test_cases:
            cut = self.CUT.copy()
            isa_qc = self.pass_manager.run(cut)
            sp = [(k, v) for k, v in test_case["test_case"].items()]
            M1 = SparsePauliOp.from_list(sp)
            isa_observables = M1.apply_layout(isa_qc.layout)
            pubs.append((isa_qc,isa_observables))

        with Session(backend=self.aer_sim) as session:
            estimator = Estimator(mode=session,options=EstimatorOptions(default_shots=self.shots))
            results = estimator.run(pubs).result()

        return [x.data.evs for x in results]


    def run_randomsearch(self,epoch=100):
        best_fit = 0
        all_fit = []
        test = {}
        for i in range(epoch):
            testcase = self.random_test_case_Z()
            test_case = self.get_test_case_theoretics_Z(testcase)
            exp = self.get_theoretical_exp_from_test_case_M3(test_case)
            obs = self.execute_test_cases([test_case])[0]
            if abs(exp - obs)>best_fit:
                best_fit = abs(exp - obs)
                test = test_case["test_case"]
            all_fit.append(abs(exp - obs))

        return best_fit,test, all_fit

    def run_mealga32(self,fam_idx,epoch,pop_size):
        max_family_size = 32

        bounds = [IntegerVar(lb=(1,)* max_family_size, ub=((2**self.CUT.num_qubits)-1,)* max_family_size, name="paulies"),
                  FloatVar(lb=(-1.,) * max_family_size, ub=(1.,) * max_family_size, name="delta")]
        

        problem = MealGAProblem(bounds=bounds,tester=self,fam_idx=fam_idx)

        model = GA.BaseGA(epoch=epoch, pop_size=pop_size)
        g_best = model.solve(problem)


        decoded = problem.decode_solution(g_best.solution)
        indexes, weights = decoded['paulies'], decoded['delta']
        paulies = get_Z_family_values(self.CUT.num_qubits, indexes)
        pauli_dict = {}
        for pauli, prob in zip(paulies, weights):
            pauli_dict[pauli] = prob
        return g_best.target.fitness, pauli_dict, model.history

    def run_mealoneplusone(self,fam_idx,epoch):
        max_family_size = 32

        bounds = [IntegerVar(lb=(1,) * max_family_size, ub=((2 ** self.CUT.num_qubits) - 1,) * max_family_size,
                             name="paulies"),
                  FloatVar(lb=(-1.,) * max_family_size, ub=(1.,) * max_family_size, name="delta")]
        

        problem = MealGAProblem(bounds=bounds, tester=self, fam_idx=fam_idx)

        model = ES.OriginalES(epoch=epoch,pop_size=5)
        g_best = model.solve(problem)

        

        decoded = problem.decode_solution(g_best.solution)
        indexes, weights = decoded['paulies'], decoded['delta']
        paulies = get_Z_family_values(self.CUT.num_qubits, indexes)
        pauli_dict = {}
        for pauli, prob in zip(paulies, weights):
            pauli_dict[pauli] = prob
        return g_best.target.fitness, pauli_dict, model.history
    
    def run_reverse_mealoneplusone(self,fam_idx,epoch):
        max_family_size = 32

        bounds = [IntegerVar(lb=(1,) * max_family_size, ub=((2 ** self.CUT.num_qubits) - 1,) * max_family_size,
                             name="paulies"),
                  FloatVar(lb=(-1.,) * max_family_size, ub=(1.,) * max_family_size, name="delta")]
        

        problem = ReverseMealGAProblem(bounds=bounds, tester=self, fam_idx=fam_idx)

        model = ES.OriginalES(epoch=epoch,pop_size=5)
        g_best = model.solve(problem)

        

        decoded = problem.decode_solution(g_best.solution)
        indexes, weights = decoded['paulies'], decoded['delta']
        paulies = get_Z_family_values(self.CUT.num_qubits, indexes)
        pauli_dict = {}
        for pauli, prob in zip(paulies, weights):
            pauli_dict[pauli] = prob
        return g_best.target.fitness, pauli_dict, model.history

    def run_mealhillclimbing(self,fam_idx,epoch,pop_size):
        max_family_size = 32

        bounds = [IntegerVar(lb=(1,) * max_family_size, ub=((2 ** self.CUT.num_qubits) - 1,) * max_family_size,
                             name="paulies"),
                  FloatVar(lb=(-1.,) * max_family_size, ub=(1.,) * max_family_size, name="delta")]
        

        problem = MealGAProblem(bounds=bounds, tester=self, fam_idx=fam_idx)

        model = HC.OriginalHC(epoch=epoch, pop_size=pop_size)
        g_best = model.solve(problem)

        decoded = problem.decode_solution(g_best.solution)
        indexes, weights = decoded['paulies'], decoded['delta']
        paulies = get_Z_family_values(self.CUT.num_qubits, indexes)
        pauli_dict = {}
        for pauli, prob in zip(paulies, weights):
            pauli_dict[pauli] = prob
        return g_best.target.fitness, pauli_dict, model.history


class MealGAProblem(MealProblem):
    def __init__(self, bounds=None, minmax="max", tester=None,fam_idx=0, **kwargs):
        self.tester = tester
        self.fam_idx = fam_idx
        super().__init__(bounds, minmax, **kwargs)

    def obj_func(self, x):
        decoded = self.decode_solution(x)
        indexes, weights = decoded['paulies'], decoded['delta']
        fam = self.tester.applicable_families[self.fam_idx]
        paulies = get_Z_family_values(self.tester.CUT.num_qubits,indexes)
        pauli_dict = {}
        for pauli, prob in zip(paulies, weights):
            pauli_dict[pauli] = prob
        testcase = {"test_case": pauli_dict, "family_index": fam[0], "M": fam[-1]}
        if testcase["test_case"] == {}:
            return np.inf
        else:
            test_case = self.tester.get_test_case_theoretics_Z(testcase)
            exp = self.tester.get_theoretical_exp_from_test_case_M3(test_case)
            obs = self.tester.execute_test_cases([test_case])[0]
        return abs(exp - obs)
    
class ReverseMealGAProblem(MealProblem):
    def __init__(self, bounds=None, minmax="min", tester=None,fam_idx=0, **kwargs):
        self.tester = tester
        self.fam_idx = fam_idx
        super().__init__(bounds, minmax, **kwargs)

    def obj_func(self, x):
        decoded = self.decode_solution(x)
        indexes, weights = decoded['paulies'], decoded['delta']
        fam = self.tester.applicable_families[self.fam_idx]
        paulies = get_Z_family_values(self.tester.CUT.num_qubits,indexes)
        pauli_dict = {}
        for pauli, prob in zip(paulies, weights):
            pauli_dict[pauli] = prob
        testcase = {"test_case": pauli_dict, "family_index": fam[0], "M": fam[-1]}
        if testcase["test_case"] == {}:
            return np.inf
        else:
            test_case = self.tester.get_test_case_theoretics_Z(testcase)
            exp = self.tester.get_theoretical_exp_from_test_case_M3(test_case)
            obs = self.tester.execute_test_cases([test_case])[0]
        return abs(exp - obs)

class Circuit_Tester_IBM_ZNE:
    def __init__(self,CUT:QuantumCircuit,simulator_type="Noise",device_name="ibm_brisbane",shots=10000):
        """
        :param CUT: A qiskit circuit without measurements
        """

        self.compact_program_specification = None
        self.CUT = CUT
        self.device_name = device_name
        service = QiskitRuntimeService(channel="ibm_cloud",token="add your own tokken")
        self.real_backend = service.backend(self.device_name)
        if simulator_type=="Noise":
            self.aer_sim = AerSimulator.from_backend(self.real_backend,method="statevector", device='GPU', blocking_enable=True, batched_shots_gpu=True)
        elif simulator_type=="Real":
            pass
        self.shots = shots
        self.simulator_type = simulator_type

    def set_applicable_families_Z(self,compact_program_specification:dict):
        """
        :param compact_program_specification: dict of paulistring with their measurement results
        """
        self.compact_program_specification = compact_program_specification
        self.applicable_families = []
        zfam = "Z"*(2**self.CUT.num_qubits)
        for key,value in compact_program_specification.items():
            if key in zfam:
                self.applicable_families.append((0, value))

    def random_test_case_Z(self):
        pauli_dict = {}
        fam = random.choice(self.applicable_families)
        selected_pauli = get_random_Z_family(self.CUT.num_qubits)
        for s in selected_pauli:
            pauli_dict[s] = np.random.random()
        return {"test_case":pauli_dict,"family_index":fam[0],"M":fam[-1]}


    def get_test_case_theoretics_Z(self,test_case:dict):
        return {"test_case": test_case["test_case"], "M": test_case["M"]}

    def get_theoretical_exp_from_testcase(self, test_case:dict):

        def pauli_expectation_from_counts(counts, pauli_string):
            pauli_string = pauli_string.upper()
            total_counts = sum(counts.values())
            exp_val = 0.0
            for bitstring, count in counts.items():
                v = 1
                for i, p in enumerate(pauli_string):
                    if p == 'Z':
                        v *= 1 if bitstring[i] == '0' else -1
                exp_val += v * count
            return exp_val / total_counts

        # Compute expectation value of the Hamiltonian
        counts = test_case["M"]
        expectation = sum(
            coeff * pauli_expectation_from_counts(counts, pauli_str)
            for pauli_str, coeff in test_case['test_case'].items()
        )
        return expectation

    def get_theoretical_exp_from_test_case_M3(self, test_case:dict):

        counts = test_case["M"]
        quasi = vector_to_quasiprobs(counts_to_vector(counts),counts)
        expectation = sum([v*quasi.expval(exp_ops=x) for x,v in test_case["test_case"].items()])
        return expectation

    def execute_test_cases(self,test_cases:list):
        from mitiq import zne
        pass_manager = generate_preset_pass_manager(backend=self.aer_sim, optimization_level=0)
        exps = []
        for test_case in test_cases:
            cut = self.CUT.copy()

            def execute_mitiq(circuit):

                isa_qc = pass_manager.run(circuit)
                sp = [(k, v) for k, v in test_case.items()]
                M1 = SparsePauliOp.from_list(sp)
                isa_observables = M1.apply_layout(isa_qc.layout)
                pubs = (isa_qc,isa_observables)
                with Session(backend=self.aer_sim) as session:
                    estimator = Estimator(mode=session,options=EstimatorOptions(default_shots=self.shots))
                    results = estimator.run([pubs]).result()
                    return [x.data.evs for x in results][0]

            def medoid(values):
                values = np.asarray(values, dtype=float)
                n = len(values)
                total_dists = [np.sum(np.abs(values[i] - values)) for i in range(n)]
                idx = int(np.argmin(total_dists))
                return values[idx]

            poly_factory = zne.inference.PolyFactory(scale_factors=[1.0, 2.0, 2.5, 3.0],order=3)
            exp = medoid([zne.execute_with_zne(cut,execute_mitiq,factory=poly_factory,scale_noise=zne.scaling.fold_global) for x in range(5)])
            exps.append(exp)
    def run_customoneplusone(self,fam_idx,epoch):
        max_family_size = 32
        space = []
        for _ in range(max_family_size//2):
            space.append({"lo": 1.0, "hi": (2 ** self.CUT.num_qubits)-1, "type": "int"})

        for _ in range(max_family_size//2):
            space.append({"lo": -1.0, "hi": 1.0, "type": "float"})

        search_space: SearchSpace = space

        es = customES(
            search_space=search_space,
            pop_size=5,
            num_children=3,
            max_gens=epoch,
            target_fitness=2,
            seed=None,
        )

        es.tester = self
        es.fam_idx = fam_idx

        best,history = es.run()
        decoded = es.decode_vector(best["vector"])
        best_fit = best['fitness']
        indexes, weights = decoded[0:len(decoded)//2], decoded[len(decoded)//2::]
        paulies = get_Z_family_values(self.CUT.num_qubits, indexes)
        pauli_dict = {}
        for pauli, prob in zip(paulies, weights):
            pauli_dict[pauli] = prob

        return best_fit, pauli_dict, history



if __name__ == '__main__':
    pass
