# Created by roy.gonzalez-aleman at 09/03/2024
"""
Plot general insights from the benchmark
"""
import os
from os.path import join

import pandas as pd
import prody as prd

import commons as cmn


def get_last_n_lines(file_name, N):
    """
    Emulate the behaviour of tail command in bash

    Args:
        file_name: path to the file
        N: number of lines to read
    """
    list_of_lines = []
    with open(file_name, 'rb') as read_obj:
        read_obj.seek(0, os.SEEK_END)
        buffer = bytearray()
        pointer_location = read_obj.tell()
        while pointer_location >= 0:
            read_obj.seek(pointer_location)
            pointer_location -= 1
            new_byte = read_obj.read(1)
            if new_byte == b'\n':
                list_of_lines.append(buffer.decode()[::-1])
                if len(list_of_lines) == N:
                    return list(reversed(list_of_lines))
                buffer = bytearray()
            else:
                buffer.extend(new_byte)
        if len(buffer) > 0:
            list_of_lines.append(buffer.decode()[::-1])
    return list(reversed(list_of_lines))


def get_n_poses(poses_path):
    poses_lines = get_last_n_lines(poses_path, 500)
    for line in poses_lines[::-1]:
        if line.startswith('MODEL'):
            n_poses = int(line.split()[-1])
            return n_poses
    return None


# ==== User-defined parameters ================================================

# Benchmark data
bench_data_path = 'scripts/04_complexes_analyses/pickles/data_dict.pick'

# Clustering data
clusters_dir = './scripts/04_complexes_analyses/03_clusterings'

# Coverage data
coverage_dir = 'scripts/04_complexes_analyses/01_rec_coverage'

# Output data
out_dir = 'scripts/04_complexes_analyses/07_insights'
os.makedirs(out_dir, exist_ok=True)

# =============================================================================

bench_data = cmn.unpickle_from_file(bench_data_path)
cases = sorted(bench_data.keys())
programs = ['plants', 'autodock4', 'vina', 'smina', 'gnina', 'qvinaw']
levels = ['low', 'medium', 'high']

insights = cmn.recursive_defaultdict()
for case in cases:
    mapping = next(
        cmn.recursive_finder('mappings.txt', join(coverage_dir, case)))
    mapping_df = pd.read_table(mapping, sep='\s+')
    total_res_nums, total_num_sf = mapping_df.shape
    coverage_df = dict((mapping_df > 0).sum() / total_res_nums * 100)

    for program in programs:
        for sf in sorted(bench_data[case][program].keys()):
            for level in levels:
                print(f'Analyzing {case}-{program}-{sf}-{level}')
                # Info from 00_get_all_data
                info = bench_data[case][program][sf][level]
                poses = info['poses_filtered']
                ag = prd.Ensemble(prd.parsePDB(poses))
                n_poses = get_n_poses(poses)

                # Info from coverage
                coverage = coverage_df[f'{program}_{sf}_{level}']

                # Perform clustering
                n_int = len(cmn.leader_clustering_traj(poses, None, 0.2))
                n_ext = len(cmn.leader_clustering_matrix(ag, 2))

                # Populate insights
                insights[case][program][sf][level]['n_poses'] = n_poses
                insights[case][program][sf][level]['n_clust_int'] = n_int
                insights[case][program][sf][level]['n_clust_ext'] = n_ext
                insights[case][program][sf][level]['coverage'] = coverage
                insights[case][program][sf][level]['time'] = info['time']
                insights[case][program][sf][level]['ram'] = info['ram']

    cmn.pickle_to_file(insights, join(out_dir, 'insights.pickle'))
