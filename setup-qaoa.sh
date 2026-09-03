#!/bin/bash -l

#create an environment

conda create --name qaoa-setup python=3.9.6
conda init bash
conda activate qaoa-setup

#package installation

python3.9 -m pip install -U qiskit==0.45.2
python3.9 -m pip install -U qiskit-aer==0.15.1
python3.9 -m pip install -U qiskit-algorithms==0.3.1
python3.9 -m pip install -U qiskit-optimization==0.6.1

python3.9 -m pip install openqaoa-braket==0.2.6
python3.9 -m pip install openqaoa-core==0.2.6
python3.9 -m pip install openqaoa-qiskit==0.2.6


python3.9 -m pip install cachetools==5.5.0
python3.9 -m pip install pennylane==0.38.0
