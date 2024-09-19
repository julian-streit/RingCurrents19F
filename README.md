# RingCurrents19F
This python script calculates distances and angles between fluorine atoms and aromatic rings to predict suitable labelling sites for 19F NMR experiments.

# Requirements
Python3.7 or later

# Installation
```
git clone https://github.com/julian-streit/RingCurrents19F.git
cd RingCurrents19F
python -m venv pyRingCurrents19F
source pyRingCurrents19F/bin/activate
pip install numpy MDAnalysis
python ring_19F.py -h
```

# Testing/Example
```
mkdir output_dir
mkdir AF_output_dir
```
Example of running the script on an MD trajectory:
```
