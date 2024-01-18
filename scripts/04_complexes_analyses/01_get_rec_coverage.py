# Created by roy.gonzalez-aleman at 02/01/2024
import os
from os.path import join

import pandas as pd
from scipy.spatial import cKDTree as ckd

import commons as cmn
import root

# todo: smooth color for visual please
# todo: put a row with normalized distance from the active site
# todo: recover the active site residue in the preparation step


# ==== User-defined parameters ================================================
all_data_pick = 'scripts/04_complexes_analyses/pickles/data_dict.pick'
inputs_dir = 'scripts/03_complexes_explorations/input_files'
top_out_dir = 'scripts/04_complexes_analyses/01_rec_coverage'
os.makedirs(top_out_dir, exist_ok=True)
# =============================================================================


# ==== Synonyms dictionaries

# ==== Orderings
programs_order = ['plants', 'qvinaw', 'vina', 'smina']

# ==== Start processing
all_data = cmn.unpickle_from_file(all_data_pick)
for case in all_data:
    out_dir = join(top_out_dir, case)
    os.makedirs(out_dir, exist_ok=True)
    # Process a single receptor conformation
    input_dir = join(inputs_dir, case)
    rec_path = next(cmn.recursive_finder('receptor.pdb', input_dir))
    rec_parsed = root.Molecule(rec_path).parse()[0]
    rec_res_nums = rec_parsed.getResnums()
    rec_ckd = ckd(rec_parsed.getCoords())

    # Start iterating for mapping update
    mapping = dict()
    for program in sorted(all_data[case].keys()):
        for sf in sorted(all_data[case][program].keys()):
            for level in ['low', 'medium', 'high']:
                poses_filtered = all_data[case][program][sf][level][
                    'poses_filtered']

                # Detect rec-lig contacts under r A
                mapping_name = f'{program}_{sf}_{level}'
                mapping.update({mapping_name: {}})
                for i in range(poses_filtered.numCoordsets()):
                    poses_filtered.setACSIndex(i)
                    pose_ckd = ckd(poses_filtered.hydrogen.getCoords())
                    contacts_nested = pose_ckd.query_ball_tree(rec_ckd, r=5)
                    contacts_atoms = [y for x in contacts_nested for y in x]
                    contacts_res_nums = set(rec_res_nums[contacts_atoms])

                    # Update mapping
                    for res_num in contacts_res_nums:
                        if res_num in mapping[mapping_name]:
                            mapping[mapping_name][res_num] += 1
                        else:
                            mapping[mapping_name].update({res_num: 1})
    df = pd.DataFrame(mapping)
    with open(join(out_dir, 'mappings.txt'), 'wt') as txt:
        df.to_string(txt)
