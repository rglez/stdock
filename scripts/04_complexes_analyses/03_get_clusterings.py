# Created by roy.gonzalez-aleman at 02/01/2024
import os
from os.path import join

import mdtraj as md
import numpy as np
import prody as prd
from numba import jit, prange

import commons as cmn

# =============================================================================
# User-specified parameters
# =============================================================================
all_data_pick = 'scripts/04_complexes_analyses/pickles/data_dict.pick'
all_data = cmn.unpickle_from_file(all_data_pick)
top_out_dir = 'scripts/04_complexes_analyses/03_clusterings'
os.makedirs(top_out_dir, exist_ok=True)
rmsd_cutoff1 = 0.1  # in nm, so 1 A
rmsd_cutoff2 = 2  # in A


@jit(nopython=True, parallel=True)
def rmsd_ref_vs_all(ref, array, divByN):
    rmsds = np.zeros(len(array), dtype=float)
    for i in prange(len(array)):
        tar = array[i]
        rmsds[i] = np.sqrt(((ref - tar) ** 2).sum() * divByN)
    return rmsds


def reorder_poses(poses, ordering):
    new_poses = prd.AtomGroup()
    for pose in poses.iterCoordsets():
        new_poses.addCoordset(pose[ordering])
    return new_poses


def leader_clustering_matrix(sorted_ensemble, cutoff):
    clusters = []
    N = len(sorted_ensemble)
    clustered = np.zeros(N, dtype=bool)
    array = sorted_ensemble.getCoordsets()
    divByN = 1.0 / array[0].shape[0]

    while not clustered.all():
        next_frame = clustered.argmin()
        rmsd = rmsd_ref_vs_all(array[next_frame], array, divByN)
        new_cluster = rmsd <= cutoff
        true_clusters = np.bitwise_and(new_cluster, ~clustered)
        clustered[new_cluster] = True
        clusters.append(true_clusters.nonzero()[0])
    return clusters


def leader_clustering_traj(sorted_pdb_path, cutoff):
    trajectory = md.load(sorted_pdb_path)
    trajectory.center_coordinates()

    clusters = []
    N = trajectory.n_frames
    clustered = np.zeros(N, dtype=bool)

    while not clustered.all():
        next_frame = clustered.argmin()
        rmsd = md.rmsd(trajectory, trajectory, next_frame)
        new_cluster = rmsd <= cutoff
        true_clusters = np.bitwise_and(new_cluster, ~clustered)
        clustered[new_cluster] = True
        clusters.append(true_clusters.nonzero()[0])
    return clusters


for case in all_data:
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
                # Extract poses
                poses = all_data[case][prog][sf][level]['poses_filtered']
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
    sorted_ensemble_path = prd.writePDB(
        join(out_dir, f'sorted_poses.pdb'), sorted_ensemble)

    # Clustering internal (superposing the ensemble)
    clusters_super = leader_clustering_traj(sorted_ensemble_path, rmsd_cutoff1)
    cmn.pickle_to_file(clusters_super, join(out_dir, 'clusters_super.pick'))

    # Clustering external (without superposing the ensemble)
    sorted_ensemble.setCoords(sorted_ensemble[0])
    clusters_nosuper = leader_clustering_matrix(sorted_ensemble, rmsd_cutoff2)
    cmn.pickle_to_file(clusters_nosuper,
                       join(out_dir, 'clusters_nosuper.pick'))
