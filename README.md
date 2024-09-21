# RingCurrents19F
This python script calculates distances and angles between fluorine atoms and aromatic rings to predict suitable labelling sites for 19F NMR experiments.

# Requirements
Python3.10 or later

# Installation
```
git clone https://github.com/julian-streit/RingCurrents19F.git
cd RingCurrents19F
python -m venv pyRingCurrents19F
source pyRingCurrents19F/bin/activate
pip install numpy MDAnalysis
python ring_19F.py -h
deactivate
```

# Testing/Examples
```
mkdir output_dir
mkdir AF_output_dir
source pyRingCurrents19F/bin/activate
```
Example of running the script on an MD trajectory:
```
python ring_19F.py -pdb MD_data/topology.pdb -xtc MD_data/traj_short.xtc -stride 1 -res1 655 -res2 675 -atoms CH FH1 FH2 FH3 -output output_dir/ -geomfactor True
```

Example of running the script on a single pdb file:
```
python ring_19F.py -pdb MD_data/topology.pdb -res1 655 -res2 675 -atoms CH FH1 FH2 FH3 -output output_dir/ -geomfactor True
```

Example of running the script on multiple pdb files stored in a directory (batch analysis):
```
python ring_19F.py -dir AF_data/ -res1 11 -res2 31 -atoms OH -output AF_output_dir/ -geomfactor True
```
Deactivate the virtual environment when finished:
```
deactivate
```
