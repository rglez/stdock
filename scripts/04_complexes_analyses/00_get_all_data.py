# Created by roy.gonzalez-aleman at 02/01/2024
import os
import re
from os.path import join

import numpy as np
import pandas as pd

import commons as cmn


# todo: ensure correct ordering will be used in rmsd calculations (AD4 issue)


def get_time_and_ram(log_file_path):
    string_time = 'Elapsed.*:*'
    string_ram = 'Maximum resident.*:*'
    with open(log_file_path, 'rt') as log:
        file_string = log.read()
    time_raw = re.findall(string_time, file_string)[0]
    time = time_raw.split(': ')[-1].split(':')
    if len(time) == 2:
        real_time = int(time[0]) * 60 + float(time[1])
    elif len(time) == 3:
        real_time = int(time[0]) * 3600 + int(time[1]) * 60 + float(time[2])
    else:
        raise ValueError('Problems parsing Elapsed time')

    ram_raw = re.findall(string_ram, file_string)[0]
    ram = int(ram_raw.split('):')[-1])
    real_ram = ram / 2 ** 20
    return real_time, real_ram


def get_scores(scores_txt_path):
    df = pd.read_table(scores_txt_path, sep='\s+', names=['frame', 'score'])
    return np.asarray(df.score)


# =============================================================================
# User-defined parameters
# =============================================================================
# exploration_dir = 'scripts/03_complexes_explorations/explorations'
exploration_dir = '/media/roy.gonzalez-aleman/Expansion/01_explored_completed/'
output_dir = 'scripts/04_complexes_analyses/pickles'

# =============================================================================
# Get all data and process it in the adequate format
# =============================================================================
data_dict = cmn.recursive_defaultdict()

cases = os.listdir(exploration_dir)
f = 0
t = 0
for case in cases:
    programs_names = os.listdir(join(exploration_dir, case))
    for program_name in programs_names:
        explorations = os.listdir(join(exploration_dir, case, program_name))

        for exploration in explorations:
            # Get scoring function & exhaustiveness
            splitted = exploration.split('_')
            if len(splitted) == 2:
                sf, level = splitted
            elif len(splitted) == 3:
                _, sf, level = splitted
            elif len(splitted) == 4:
                _, sf1, sf2, level = splitted
                sf = '_'.join([sf1, sf2])
            else:
                raise ValueError('Problems with parsing SF')

            # Save top directory
            top_dir = join(exploration_dir, case, program_name, exploration)
            data_dict[case][program_name][sf][level]['top_dir'] = top_dir

            # Save log_file information
            log_file = next(
                cmn.recursive_finder(f'{exploration}.log', top_dir))
            seconds, gb = get_time_and_ram(log_file)
            data_dict[case][program_name][sf][level]['time'] = seconds
            data_dict[case][program_name][sf][level]['ram'] = gb

            # Save all poses
            poses = next(cmn.recursive_finder('poses.pdb', top_dir))
            data_dict[case][program_name][sf][level]['poses'] = poses

            # Save filtered poses
            poses_filtered = next(
                cmn.recursive_finder('poses_filtered.pdb', top_dir))
            data_dict[case][program_name][sf][level][
                'poses_filtered'] = poses_filtered

            # Save ordering
            npy = list(cmn.recursive_finder('*.npy', top_dir))
            if npy:
                data_dict[case][program_name][sf][level]['order'] = npy[0]

            # Save all scores
            scores = get_scores(
                next(cmn.recursive_finder('scores.txt', top_dir)))
            data_dict[case][program_name][sf][level]['scores'] = scores

            # Save filtered scores
            scores_filtered = get_scores(
                next(cmn.recursive_finder('scores_filtered.txt', top_dir)))
            data_dict[case][program_name][sf][level][
                'scores_filtered'] = scores_filtered

os.makedirs(output_dir, exist_ok=True)
cmn.pickle_to_file(data_dict, join(output_dir, 'data_dict.pick'))
