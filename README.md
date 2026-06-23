# Heisenberg SBFL Experiment
This repository serves as the code base for the UCCS QAS Lab's experiment with SBFL in the Heisenberg picture. This repository was forked and built upon the SB-QOPS framework, which can be found [here](https://github.com/AsmarMuqeet/SB-QOPS).

# Placeholder Activity Diagram
![Placeholder Activity Diagram](SBFL_assets/Placeholder%20Activity%20Diagram.png)

# Installation
The repository is only tested in a Linux environment since Qiskit AER GPU is only supported in linux

### Dependencies

Anaconda Python distribution is required [here](https://www.anaconda.com/products/distribution):

Steps:

    1. Clone the repository
    2. cd SB-QOPS
    3  conda env create -f environment.yml
    4. conda activate sbqops 
	
# Evaluate new Circuit:
### To test new circuits

``` python

import QOPS as qops

if __name__ == '__main__':

    QUBITS = 29 # number of qubits

    ga_result = pd.DataFrame(columns=['Program',"mutant",'catch_avg','avg_fitness','testcases'])
    circuit = # qiskit circuit without measurements
    program_specification = #compact program specification in the form {"paulistring eg zzzz": {"bitstrings":count}}
    tester = qops.Circuit_Tester(CUT=circuit)
    tester.set_applicable_families_Z(program_specification)
    for i in range(len(tester.applicable_families)):
        best_function,testcase, history = tester.run_mealoneplusone(i, 80) # change algorithm here
            if best_function > 0.1: # tolerance threshold
                killed = 1
                pauli = testcase
                fitness = best_function
                break
```

# Acknowledgements
We acknowledge and thank the authors of SB-QOPS Asmar Muqeet, Shaukat Ali, and Paolo Arcaini for their work and making their work open-source under the Apache 2.0 License