# Created by roy.gonzalez-aleman at 24/12/2023
"""
This script is conceived to cluster the X conformers generated for each ligand
 in the DUDE dataset and then select the most mobile as case study.
"""
import os
from os.path import join, split

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

import commons as cmn


def get_n_rot(input_pdb):
    mol = Chem.AddHs(Chem.MolFromPDBFile(input_pdb), addCoords=True)
    return Chem.rdMolDescriptors.CalcNumRotatableBonds(mol), mol.GetAtoms()


# =============================================================================
# User-specified parameters
# =============================================================================
root_dir = os.path.abspath('.')
conformers_dir = join(root_dir,
                      'scripts/dude_selection/01_conformers_generation')
output_dir = join(root_dir, 'scripts/dude_selection/02_clustering')
os.makedirs(output_dir, exist_ok=True)
min_clust_size = 2
cutoff = 2
bitqt = '/home/roy.gonzalez-aleman/miniconda3/envs/stdock/bin/bitqt'

# =============================================================================
# Clustering
# =============================================================================

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

cmn.pickle_to_file(topo_trajs_dict, join(root_dir, 'scripts/dude_selection/topo_traj.pick'))

# Run clustering jobs with mdscan
for case in topo_trajs_dict:
    topo = topo_trajs_dict[case]['top']
    traj = topo_trajs_dict[case]['traj']
    out_name = join(output_dir, f'{case}_bitqt')
    os.system(
        f'{bitqt} -top {topo} -traj {traj} -odir {out_name} -cutoff {cutoff} -min_clust_size {min_clust_size}')
