#!/bin/python

# Import modules
import numpy as np
import MDAnalysis.analysis.distances
import MDAnalysis as md
import argparse
import os

# define command line input
parser = argparse.ArgumentParser(prog='Ring-19F',description='This python script calculates distances and angles between fluorine atoms and aromatic rings',
    epilog='Units: Distances = Angstrom, Angles = Degrees, Geomectric factors = cm^-3')

parser.add_argument('-pdb',default=None,help='PDB file of structure to be analysed and/or topology for MD trajectory')
parser.add_argument('-dir',default=None,help='directory containing multiple pdb files for batch analysis')
parser.add_argument('-xtc',default=None,help='XTC trajectory file (has to be compatible with the PDB option)')
parser.add_argument('-stride',default=1,type=int,help='Stride for trajectory analysis (e.g., stride=10 means that every 10th frame is analysed)')
parser.add_argument('-res1',required=True,help='19F labelled residue number')
parser.add_argument('-res2',required=True,help='Aromatic residue number')
parser.add_argument('-atoms',required=True,nargs='+',help='Fluorine atom name(s) - if multiple are provided the center of mass of the listed atoms will be used')
parser.add_argument('-output',required=True,help='Output directory')
parser.add_argument('-geomfactor',default=False,help='Set to True to get the geometric factor, (1-3cos^2(theta) / r^3')
args = parser.parse_args()

# define functions
def unit_vector(vector): # returns the unit vector of a vector

    return vector / np.linalg.norm(vector)

def angle_between(v1, v2): # angle between vectors v1 and v2

    v1_u = v1
    v2_u = v2
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def calc_distance(u,res_19F,aromatic_res,aromatic_type,atom_names_19F,stride=1):

    # initiate list for distances
    d = []
    frames = np.arange(len(u.trajectory))[::stride] # frames of trajectory (for pdb files this is just 1 frame)

    # atom selections
    atoms19F = str('name {}'.format(atom_names_19F[0]))
    if len(atom_names_19F)>1:
        for at in atom_names_19F[1:]:
            atoms19F+=str(' or name {}'.format(at))
    sel_19F = str('resid {} and ({})'.format(res_19F,atoms19F))

    ring_sel = str('resid {} and not backbone and not name H* and not name CB and not name OH'.format(aromatic_res))

    # iterate through all frames and calculate distances
    for i in frames:
        u.trajectory[i]
        c1 = u.select_atoms(sel_19F).center_of_mass()
        c2 = u.select_atoms(ring_sel).center_of_mass()
        dist = MDAnalysis.analysis.distances.distance_array(c1,c2)[0][0]
        d.append(dist)

    return np.array(d)

def calc_angle(universe,res_19F,aromatic_res,aromatic_type,atom_names_19F,stride=1):

    # initiate list for distances
    angles = []
    frames = np.arange(len(u.trajectory))[::stride] # frames of trajectory (for pdb files this is just 1 frame)

    # atom selections for calculating the ring normal vector
    sel_dict = {}
    sel_dict['PHE'] = 'and (name CG or name CE1 or name CE2)'
    sel_dict['TYR'] = 'and (name CG or name CE1 or name CE2)'
    sel_dict['HIS'] = 'and (name CG or name NE2 or name ND1)'
    sel_dict['HIE'] = 'and (name CG or name NE2 or name ND1)'
    sel_dict['HID'] = 'and (name CG or name NE2 or name ND1)'
    sel_dict['TRP'] = 'and (name CE2 or name CZ3 or name CG)'
    plane_sel = sel_dict[aromatic_type]

    # ring center of mass selection
    com_sel = 'and not backbone and not name H* and not name CB and not name OH'

    # atom selections
    atoms19F = str('name {}'.format(atom_names_19F[0]))
    if len(atom_names_19F)>1:
        for at in atom_names_19F[1:]:
            atoms19F+=str(' or name {}'.format(at))
    sel_19F = str('resid {} and ({})'.format(res_19F,atoms19F))

    # iterate through all frames and analyse angles
    for i in frames:
        u.trajectory[i]
        c1 = u.select_atoms(sel_19F).center_of_mass() # CF3 group

        # sidechain ring
        s2 = u.select_atoms('resid {} '.format(aromatic_res)+plane_sel).positions # shape Natoms x 3
        com = u.select_atoms('resid {} '.format(aromatic_res)+com_sel).center_of_mass() # shape 1x3

        # now calculate vectors p1 to p2 and p1 to p3 - these two vectors define the plane of the ring
        v1 = s2[1,:]-s2[0,:]
        v2 = s2[2,:]-s2[0,:]
        v1 = unit_vector(v1)
        v2 = unit_vector(v2)

        # now calculate a normal vector - perpendicular to the plane (i.e. cross product of v1 and v2)
        norm = unit_vector(np.cross(v1,v2))

        # now get vector of CF3 COM to ring COM
        a1 = com-c1
        a1 = unit_vector(a1)

        # now calculate angle between normal and bond vector defining the amide pi interactions
        ang1 = np.rad2deg(angle_between(a1,norm))
        angles.append(ang1)

    return np.array(angles)


# load structures / trajectories

# run analysis
# define whether this is an MD traj analysis or PDB batch run
if args.dir==None:
    if args.xtc==None:
        u = md.Universe(args.pdb)
        aromatic_residue = u.select_atoms('resid {}'.format(args.res2)).resnames[0]
    else:
        u = md.Universe(args.pdb,args.xtc)
        aromatic_residue = u.select_atoms('resid {}'.format(args.res2)).resnames[0]

    distances = calc_distance(u,args.res1,args.res2,aromatic_residue,list(args.atoms),stride=args.stride)
    angles = calc_angle(u,args.res1,args.res2,aromatic_residue,list(args.atoms),stride=args.stride)

    # save data
    np.savetxt(args.output+'/distances_{}F_{}{}.txt'.format(args.res1,args.res2,aromatic_residue),distances)
    np.savetxt(args.output+'/angles_{}F_{}{}.txt'.format(args.res1,args.res2,aromatic_residue),angles)

    if args.geomfactor=='True':
        distances_cm = distances*1e-8
        angles_rad = np.deg2rad(angles)
        geom = (1-3*np.cos(angles_rad)**2)/(distances_cm**3)

        # save data
        np.savetxt(args.output+'/geomfactors_{}F_{}{}.txt'.format(args.res1,args.res2,aromatic_residue),geom)

else:
    pdb_files = []
    path=args.dir
    for root, dir, files in os.walk(path+'/'):
        for file in files:
            if '.pdb' in file:
                pdb_files.append(file)
    pdb_files.sort() # sort pdb files

    # analyse all pdb pdb_files
    for FILE in pdb_files:
        name = FILE[:-4] # extract name by rmoving .pdb ending

        # load structure
        u = md.Universe(args.dir+str('/')+FILE)
        aromatic_residue = u.select_atoms('resid {}'.format(args.res2)).resnames[0]
        distances = calc_distance(u,args.res1,args.res2,aromatic_residue,list(args.atoms),stride=args.stride)
        angles = calc_angle(u,args.res1,args.res2,aromatic_residue,list(args.atoms),stride=args.stride)

        # save data
        np.savetxt(args.output+'/{}_distances_{}F_{}{}.txt'.format(name,args.res1,args.res2,aromatic_residue),distances)
        np.savetxt(args.output+'/{}_angles_{}F_{}{}.txt'.format(name,args.res1,args.res2,aromatic_residue),angles)

        if args.geomfactor=='True':
            distances_cm = distances*1e-8
            angles_rad = np.deg2rad(angles)
            geom = (1-3*np.cos(angles_rad)**2)/(distances_cm**3)

            # save data
            np.savetxt(args.output+'/{}_geomfactors_{}F_{}{}.txt'.format(name,args.res1,args.res2,aromatic_residue),geom)


print('Done!')
