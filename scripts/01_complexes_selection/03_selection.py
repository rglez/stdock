# Created by roy.gonzalez-aleman at 24/12/2023
import os
import shutil
from collections import defaultdict
from os.path import join, basename

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

import commons as cmn

# =============================================================================
# Selection
# =============================================================================

# ==== Prepare folders hierarchy
root_dir = cmn.proj_path
topo_traj_pick = join(
    root_dir, 'scripts/01_complexes_selection/topo_traj.pick')
out_dir = join(root_dir, 'scripts/01_complexes_selection/03_selection')
clusters_dir = join(root_dir, 'scripts/01_complexes_selection/02_clustering')
shutil.rmtree(out_dir, ignore_errors=True)
os.makedirs(out_dir)

# ==== Get number of clusters produced in the previous step
topo_traj = cmn.unpickle_from_file(topo_traj_pick)
pickles = list(cmn.recursive_finder('*.pick', clusters_dir))
for pickle in pickles:
    label = basename(pickle).split('_')[1]
    clust_array = cmn.unpickle_from_file(pickle)
    n_clusters = clust_array.max()
    topo_traj[label]['num_clusters'] = n_clusters

# ==== Classify by number of rotatable bonds
by_nrots = defaultdict(list)
for case in topo_traj:
    by_nrots[topo_traj[case]['n_rot']].append(topo_traj[case])

# ==== Start selectioin procedure
selected = {}
for n_rot in by_nrots:
    options = by_nrots[n_rot]
    max_n_clusters = 0
    max_n_atoms = 0

    # for every option with the same number of rotatable bonds
    for option in options:
        n_clusters = option['num_clusters'] if option['num_clusters'] else 1
        n_atoms = option['n_atoms']

        # Ignore macrocycles
        mol = Chem.AddHs(Chem.MolFromPDBFile(option['top']), addCoords=True)
        mol_ri = mol.GetRingInfo().AtomRings()
        macrocycle = any([len(x) > 6 for x in mol_ri])
        if macrocycle:
            continue

        # Get the one with most clusters
        elif n_clusters > max_n_clusters:
            max_n_clusters = n_clusters
            max_n_atoms = n_atoms
            selected.update({n_rot: option})

        # At equal number of clusters, get the one with most atoms
        elif (n_clusters == max_n_clusters) and (n_atoms > max_n_atoms):
            max_n_clusters = n_clusters
            max_n_atoms = n_atoms
            selected.update({n_rot: option})

# =============================================================================
# Saving figs
# =============================================================================
for case in by_nrots:
    mol_suppl = [Chem.AddHs(Chem.MolFromPDBFile(x['top']), addCoords=True)
                 for x in by_nrots[case]]
    labels = [basename(x['top']).split("_")[1] for x in by_nrots[case]]
    n_clusts = [x['num_clusters'] for x in by_nrots[case]]
    n_atoms = [x['n_atoms'] for x in by_nrots[case]]
    legends = [
        f'{labels[i].upper()}    NCLUST:{n_clusts[i]}     N_ATOMS:{n_atoms[i]}'
        for i in range(len(labels))]

    for m in mol_suppl:
        tmp = AllChem.Compute2DCoords(m)
    img = Chem.Draw.MolsToGridImage(mol_suppl, molsPerRow=4,
                                    subImgSize=(400, 400), legends=legends)
    img.save(join(out_dir, f'molgrid_{case}.png'))

# =============================================================================
# Saving file
# =============================================================================
with open(join(out_dir, 'report_selection.txt'), 'wt') as rep:
    for case in sorted(selected.keys()):
        rep.write(f'{case}: {selected[case]["top"]}\n')
