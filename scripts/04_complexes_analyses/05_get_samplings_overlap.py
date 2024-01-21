# Created by roy.gonzalez-aleman at 14/01/2024
from os.path import join

import numpy as np
import pandas as pd
from numba import jit, prange
from tqdm import tqdm

import commons as cmn

# =============================================================================
# User-defined parameters
# =============================================================================
cutoff = 2


# =============================================================================


@jit(nopython=True, parallel=True)
def rmsd_ref_vs_all(ref, array, division_constant):
    """
    Computes the rmsd of reference versus all targets in array

    Args:
        ref: a numpy array of shape (N, 3)
        array: a numpy array of shape (M, N, 3)
        division_constant: a constant computed as 1 / N

    Returns:
        values: a numpy array of rmsd values
    """
    M = len(array)
    values = np.zeros(M)
    for idx in prange(M):
        rmsd_val = np.sqrt(((ref - array[idx]) ** 2).sum() * division_constant)
        values[idx] = rmsd_val
    return values


# ==== Build hierarchy
proj_dir = cmn.proj_dir
data_path = cmn.get_abs('scripts/04_complexes_analyses/pickles/data_dict.pick')
out_dir = cmn.get_abs('scripts/04_complexes_analyses/05_sampling_overlaps')
cmn.makedir_after_overwriting(out_dir)

# ==== Get samplings to compare
data_dict = cmn.unpickle_from_file(data_path)
samplings = cmn.defaultdict()
for case in data_dict:
    for program in sorted(data_dict[case].keys()):
        for sf in sorted(data_dict[case][program].keys()):
            for level in ['low', 'medium', 'high']:
                poses_filtered = data_dict[case][program][sf][level][
                    'poses_filtered']
                sampling_name = f'{program}_{sf}_{level}'
                samplings[sampling_name] = poses_filtered

# ==== Compare samplings
N = samplings[list(samplings.keys())[0]].getCoords().shape[0]
div_constant = 1.0 / N

M = len(samplings)
overlaps_matrix = np.zeros((M, M), dtype=int)
ref_counts_matrix = overlaps_matrix.copy()

# for all samplings i (ref)
for i, sampling_ref in enumerate(tqdm(samplings, total=len(samplings))):
    ref_traj = samplings[sampling_ref]
    num_frames_ref = ref_traj.numCoordsets()

    # for all samplings j (tar)
    for j, sampling_tar in enumerate(samplings):
        # if sampling i != sampling j
        if sampling_ref != sampling_tar:
            tar_traj = samplings[sampling_tar]
            tar_frames = tar_traj.getCoordsets()

            # count all times ref frames are present in tar frames
            ij_counts = 0
            for frame in range(num_frames_ref):
                ref_traj.setACSIndex(frame)
                ref_frame = ref_traj.getCoords()
                rmsds = rmsd_ref_vs_all(ref_frame, tar_frames, div_constant)
                ij_counts += (rmsds <= cutoff).any()

            # save counts to a matrix
            overlaps_matrix[i, j] = ij_counts
            ref_counts_matrix[i, j] = num_frames_ref

# ==== Save overlaps
overlaps_percent = overlaps_matrix / ref_counts_matrix * 100
pickle_name = join(out_dir,  f'overlaps_percent_at_{cutoff}A.pickle')
cmn.pickle_to_file(overlaps_percent,pickle_name)
