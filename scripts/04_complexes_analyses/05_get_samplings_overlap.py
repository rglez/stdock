# Created by roy.gonzalez-aleman at 14/01/2024

import numpy as np
import pandas as pd
from numba import jit, prange
from tqdm import tqdm

import commons as cmn


@jit(nopython=True, parallel=True)
def rmsd_ref_vs_all(ref, array, divByN, cutoff):
    values = np.zeros(len(array))
    for i in prange(len(array)):
        rmsd_val = np.sqrt(((ref - array[i]) ** 2).sum() * divByN)
        values[i] = rmsd_val
    return values <= cutoff


# =============================================================================
# User-defined parameters
# =============================================================================
data_path = '/home/roy.gonzalez-aleman/RoyHub/stdock/scripts/04_complexes_analyses/pickles/data_dict.pick'
out_dir = '/home/roy.gonzalez-aleman/RoyHub/stdock/scripts/04_complexes_analyses/05_sampling_overlaps'
cutoff = 2

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

# ==== Get
# ==== Compare samplings
N = samplings[list(samplings.keys())[0]].getCoords().shape[0]
divByN = 1.0 / N

M = len(samplings)
overlaps = np.zeros((M, M), dtype=int)
for i, sampling_ref in enumerate(tqdm(samplings, total=len(samplings))):
    ref_traj = samplings[sampling_ref]
    num_frames_ref = ref_traj.numCoordsets()

    for j, sampling_tar in enumerate(samplings):
        if sampling_ref != sampling_tar:

            tar_traj = samplings[sampling_tar]
            tar_frames = tar_traj.getCoordsets()
            ij_counts = 0
            for frame in range(num_frames_ref):
                ref_traj.setACSIndex(frame)
                ref_frame = ref_traj.getCoords()
                bit = rmsd_ref_vs_all(ref_frame, tar_frames, divByN, cutoff)
                ij_counts += bit.any()
            overlaps[i, j] = ij_counts

# %%===========================================================================
# Plotting
# =============================================================================

import matplotlib.pyplot as plt
plt.pcolormesh(pd.DataFrame(overlaps), cmap='turbo')
plt.colorbar()
plt.show()
plt.close()
