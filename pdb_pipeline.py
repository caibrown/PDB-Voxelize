import numpy as np
from io import StringIO
import sys
import matplotlib.pyplot as plt

np.set_printoptions(threshold=sys.maxsize)

def pdb_coordinates(pdb_filename: str) -> np.matrix:
    """ Given a PDB file, gets all the atoms and returns their coordinates
    as an n*3 matrix """
    with open(pdb_filename) as pdb_file:
        pdb_file_text = pdb_file.read()
    atoms = []
    for line in pdb_file_text.split('\n'):
        if line[0:4] == 'ATOM':
            x = float(line[31:38].strip(' '))
            y = float(line[39:46].strip(' '))
            z = float(line[47:54].strip(' '))
            atoms.append([x, y, z])
    rbd = atoms[318:525] # residues 319-526 = SARS-CoV-2 recepter binding domain
    return np.matrix(rbd)

def backbone(pdb_filename: StringIO, chain="all chains", CA_only=True):
    # Uses a PDB file to get the coordinates of atoms in the backbone of a protein chain
    with open(pdb_filename) as pdb_file:
        pdb_data = pdb_file.read()
    atoms = []
    backbone_atoms = ['CA'] if CA_only else ['CA', 'C', 'N']
    for line in pdb_data.split('\n'):
        # By default, look at all chains, unless just one is specified
        if line[0:4] == 'ATOM' and (chain is "all chains" or line[21] == chain) and line[13:16].rstrip(' ') in backbone_atoms:
            residue_num = int(line[24:26])
            x = float(line[31:38].strip(' '))
            y = float(line[39:46].strip(' '))
            z = float(line[47:54].strip(' '))
            atoms.append([x, y, z])
    rbd = atoms[318:525] # residues 319-526
    return np.matrix(rbd)

def centroid(atom_coords: np.matrix) -> np.matrix:
    """ Given the coordinates of each atom in a protein, as an n*3 matrix,
    returns the coordinates of the centroid as a 1*3 matrix """
    moment = sum(atom_coords)
    return moment / len(atom_coords)

def voxelize(atom_coords: np.matrix, centroid: np.matrix) -> np.matrix:
    """ Given the coordinates of each atom in a protein, as an n*3 matrix,
    and the coordinates of the centroid as a 1*3 matrix,
    returns a voxelized rendering of the protein """
    # Calculate size of voxel
    x = max(atom_coords[:,0])
    y = max(atom_coords[:,1])
    z = max(atom_coords[:,2])
    dim = int(max(x, y, z))
    voxel_cube = np.ndarray([dim, dim, dim], 'int')
    voxel_center = [int(dim/2), int(dim/2), int(dim/2)]
    for coordinate in atom_coords:
        print(int(coordinate), file=open('output.txt', 'a'))
        # Move centroid to origin,
        coordinate -= centroid
        # and align with voxel space centroid
        coordinate += voxel_center
        x = int(coordinate[0,0])
        y = int(coordinate[0,1])
        z = int(coordinate[0,2])
        # Update voxel: 1 - atom is present, 0 - no atom present
        voxel_cube[x, y, z] = 1
    return voxel_cube

# p = backbone('7a94.pdb')
# q = pdb_coordinates('7a94.pdb')
# v = voxelize(p,centroid(p))