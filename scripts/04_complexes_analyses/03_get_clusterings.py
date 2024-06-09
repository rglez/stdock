# Created by roy.gonzalez-aleman at 02/01/2024
import os
from os.path import join

import numpy as np
import prody as prd

import commons as cmn

# =============================================================================
# User-specified parameters
# =============================================================================
all_data_pick = 'scripts/04_complexes_analyses/pickles/data_dict.pick'
all_data = cmn.unpickle_from_file(all_data_pick)
top_out_dir = 'scripts/04_complexes_analyses/03_clusterings'
os.makedirs(top_out_dir, exist_ok=True)
rmsd_cutoff1 = 0.2  # in nm, so 2 A
rmsd_cutoff2 = 2  # in A


def reorder_poses(poses, ordering):
    new_poses = prd.AtomGroup()
    for pose in poses.iterCoordsets():
        new_poses.addCoordset(pose[ordering])
    return new_poses


for case in ['3udh', '3fur', '2qbp', '4gid']:
    # Create the output dir per case
    out_dir = join(top_out_dir, case)
    os.makedirs(out_dir, exist_ok=True)

    # Get the ensemble of poses & the corresponding ids array
    identifiers = []
    indices = []
    ensemble = prd.Ensemble()
    for prog in sorted(all_data[case]):
        for sf in sorted(all_data[case][prog]):
            for level in sorted(all_data[case][prog][sf]):
                print(f'Analyzing {case}:{prog}:{sf}:{level}')
                # Extract poses
                poses = prd.parsePDB(
                    all_data[case][prog][sf][level]['poses_filtered']).noh
                ordering_path = all_data[case][prog][sf][level]['order']
                if isinstance(ordering_path, str):
                    ordering = np.load(ordering_path)
                    poses = reorder_poses(poses, ordering)
                ensemble.addCoordset(poses)
                # Extract identifiers
                N = poses.numCoordsets()
                ids = [f'{prog}_{sf}_{level}_{x}' for x in range(N)]
                indices.extend(range(N))
                identifiers.extend(ids)
    identifiers = np.asarray(identifiers)
    indices = np.asarray(indices)

    # Sorting from best ranked
    sort = indices.argsort()
    sorted_identifiers = identifiers[sort]
    np.save(join(out_dir, 'sorted_identifiers'), sorted_identifiers)

    sorted_ensemble = ensemble[sort]
    sorted_ensemble.setAtoms(poses)
    sorted_ensemble_path_traj = prd.writeDCD(
        join(out_dir, f'sorted_poses.dcd'), sorted_ensemble)
    sorted_ensemble_path_top = prd.writePDB(
        join(out_dir, f'sorted_poses.pdb'), sorted_ensemble, csets=0)

    # Clustering internal (superposing the ensemble)
    clusters_super = cmn.leader_clustering_traj(sorted_ensemble_path_traj,
                                                sorted_ensemble_path_top,
                                                rmsd_cutoff1)

    cmn.pickle_to_file(clusters_super, join(out_dir, 'clusters_super.pick'))

    # Clustering external (without superposing the ensemble)
    sorted_ensemble.setCoords(sorted_ensemble[0])
    clusters_nosuper = cmn.leader_clustering_matrix(sorted_ensemble,
                                                    rmsd_cutoff2)
    cmn.pickle_to_file(clusters_nosuper,
                       join(out_dir, 'clusters_nosuper.pick'))
