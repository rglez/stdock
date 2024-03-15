# Created by roy.gonzalez-aleman at 02/01/2024
import os
from os.path import join

import pandas as pd
import prody as prd
from scipy.spatial import cKDTree as ckd

import commons as cmn

# ==== User-defined parameters ================================================
all_data_pick = 'scripts/04_complexes_analyses/pickles/data_dict.pick'
inputs_dir = 'scripts/02_complexes_preparation/01_prepared'
top_out_dir = 'scripts/04_complexes_analyses/01_rec_coverage'
os.makedirs(top_out_dir, exist_ok=True)
# =============================================================================



# ==== Start processing
all_data = cmn.unpickle_from_file(all_data_pick)
for case in all_data:
    out_dir = join(top_out_dir, case)
    os.makedirs(out_dir, exist_ok=True)
    # Process a single receptor conformation
    input_dir = join(inputs_dir, case)
    rec_path = join(inputs_dir, case, f'{case}_protein.pdb')
    rec_parsed = cmn.Molecule(rec_path).parse()[0].protein.noh
    rec_res_nums = rec_parsed.getResnums()
    rec_ckd = ckd(rec_parsed.getCoords())

    # Start iterating for mapping update
    mapping = dict()
    for program in ['plants', 'autodock4', 'vina', 'smina', 'gnina', 'qvinaw']:
        for sf in sorted(all_data[case][program].keys()):
            for level in ['low', 'medium', 'high']:
                print(f'Processing case {case}:{program}:{sf}:{level}')
                poses_filtered = prd.parsePDB(
                    all_data[case][program][sf][level]['poses_filtered']).noh

                # Detect rec-lig contacts under r A
                mapping_name = f'{program}_{sf}_{level}'
                mapping.update({mapping_name: {}})
                for i in range(poses_filtered.numCoordsets()):
                    poses_filtered.setACSIndex(i)
                    pose_ckd = ckd(poses_filtered.getCoords())
                    contacts_nested = pose_ckd.query_ball_tree(rec_ckd, r=5)
                    contacts_atoms = [y for x in contacts_nested for y in x]
                    contacts_res_nums = set(rec_res_nums[contacts_atoms])

                    # Update mapping
                    for res_num in contacts_res_nums:
                        if res_num in mapping[mapping_name]:
                            mapping[mapping_name][res_num] += 1
                        else:
                            mapping[mapping_name].update({res_num: 1})
                del poses_filtered
    df = pd.DataFrame(mapping)
    with open(join(out_dir, 'mappings.txt'), 'wt') as txt:
        df.to_string(txt)
