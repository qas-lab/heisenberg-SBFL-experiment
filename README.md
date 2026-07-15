# Heisenberg SBFL Experiment
This repository serves as the code base for the UCCS QAS Lab's experiment with SBFL in the Heisenberg picture. This repository was forked and built upon the SB-QOPS framework, which can be found [here](https://github.com/AsmarMuqeet/SB-QOPS).

# Placeholder Activity Diagram
![Placeholder Activity Diagram](SBFL_assets/Mid-Level%20HBFL%20Diagram.png)

# Installation
The repository is only tested in a Linux environment since Qiskit AER GPU is only supported in linux

### Dependencies

Anaconda Python distribution is required [here](https://www.anaconda.com/products/distribution):

Steps:

    1. Clone the repository
    2. cd heisenberg-SBFL-experiment
    3  conda env create -f environment.yml
    4. conda activate sbqops 
    5. pip install "qiskit<1.5.0"
    6. pip uninstall qiskit-aer qiskit-aer-gpu
    7. pip install --upgrade qiskit-aer
    8. pip install --upgrade qiskit-aer-gpu
    9. jupyter notebook SBFL.ipynb
	
# Evaluate new Circuit:
### To test new circuits

Please see further instructions in the SBFL.ipynb notebook file. It has been written to ease new users through the process of running our experiment.

# Acknowledgements
We acknowledge and thank the authors of SB-QOPS Asmar Muqeet, Shaukat Ali, and Paolo Arcaini for their work and making their work open-source under the Apache 2.0 License