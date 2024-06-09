# Created by roy.gonzalez-aleman at 24/12/2023
"""
Clusters the X conformers generated for each ligand in the previous step
and then select the most mobile (more clusters) as case study.

Known Issues:
    1. pandas will produce a harmless Warning on "setting with copy"
    2. bitqt raises errors for cases with a single frame
"""

import os
import shutil
from os.path import join, split

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

import commons as cmn
import proj_paths as pp


def get_n_rot(input_pdb):
    """
    Get the number of rotatable bonds

    Args:
        input_pdb: a pdb-formatted input molecule

    Returns:
        num_rot: number of rotatable bonds

    """
    mol = Chem.AddHs(Chem.MolFromPDBFile(input_pdb), addCoords=True)
    num_rot = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol), mol.GetAtoms()
    return num_rot


# =============================================================================
# User-specified parameters
# =============================================================================
rmsd_cutoff = 2
min_clust_size = 2
# =============================================================================

# ==== Prepare folders hierarchy
root_dir = pp.proj_dir
conformers_dir = join(root_dir,
                      'scripts/01_complexes_selection/01_conformers_generation')
output_dir = join(root_dir, 'scripts/01_complexes_selection/02_clustering')
shutil.rmtree(output_dir, ignore_errors=True)
os.makedirs(output_dir)
bitqt = join(root_dir, 'programs/bitqt.py')

# ==== Clustering

# Organize topologies & trajectories
topo_trajs_files = list(cmn.recursive_finder('*.pdb', conformers_dir))
topo_trajs_dict = cmn.recursive_defaultdict()

for file in topo_trajs_files:
    root, base = split(file)
    splitted = base.split('_')
    if splitted[0] == 'top':
        topo_trajs_dict[splitted[1]]['top'] = file
        num_rot, num_atoms = get_n_rot(file)
        topo_trajs_dict[splitted[1]]['n_rot'] = num_rot
        topo_trajs_dict[splitted[1]]['n_atoms'] = len(num_atoms)
    else:
        topo_trajs_dict[splitted[0]]['traj'] = file
pickle_name = join(root_dir, 'scripts/01_complexes_selection/topo_traj.pick')
cmn.pickle_to_file(topo_trajs_dict, pickle_name)

# Run clustering jobs with bitqt
for case in topo_trajs_dict:
    topo = topo_trajs_dict[case]['top']
    traj = topo_trajs_dict[case]['traj']
    out_name = join(output_dir, f'{case}_bitqt')
    os.system(
        f'{bitqt} -top {topo} -traj {traj} -odir {out_name}'
        f' -cutoff {rmsd_cutoff} -min_clust_size {min_clust_size}')
