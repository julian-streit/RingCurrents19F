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
# Usage
```
usage: Ring-19F [-h] [-pdb PDB] [-dir DIR] [-xtc XTC] [-stride STRIDE] -res1
                RES1 -res2 RES2 -atoms ATOMS [ATOMS ...] -output OUTPUT
                [-geomfactor GEOMFACTOR]

This python script calculates distances and angles between fluorine atoms and
aromatic rings

optional arguments:
  -h, --help            show this help message and exit
  -pdb PDB              PDB file of structure to be analysed and/or topology
                        for MD trajectory
  -dir DIR              directory containing multiple pdb files for batch
                        analysis
  -xtc XTC              XTC trajectory file (has to be compatible with the PDB
                        option)
  -stride STRIDE        Stride for trajectory analysis (e.g., stride=10 means
                        that every 10th frame is analysed)
  -res1 RES1            19F labelled residue number
  -res2 RES2            Aromatic residue number
  -atoms ATOMS [ATOMS ...]
                        Fluorine atom name(s) - if multiple are provided the
                        center of mass of the listed atoms will be used
  -output OUTPUT        Output directory
  -geomfactor GEOMFACTOR
                        Set to True to get the geometric factor,
                        (1-3cos^2(theta) / r^3

Units: Distances = Angstrom, Angles = Degrees, Geomectric factors = cm^-3
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
