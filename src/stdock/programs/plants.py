# Created by roy.gonzalez-aleman at 01/01/2024
import os
import subprocess as sp
from os.path import join

import pandas as pd
import prody as prd
from scipy.spatial import cKDTree as ckd
from tqdm.auto import tqdm

import commons as cmn


class Plants(cmn.Program):

    def create_config(self):
        xc, yc, zc = self.get_rec_center()
        lig_size = cmn.get_longest_component(self.lig_parsed)
        r = self.get_rec_axis().max() / 2 + lig_size

        config_paths = []
        for exhaustiveness in self.exhaustiveness_list:
            for sf in self.scoring_functions:
                suffix = f'plants_{sf}_{cmn.syn[exhaustiveness]}'
                out_dir = join(self.out_dir, suffix)
                os.makedirs(out_dir, exist_ok=True)
                config_name = join(out_dir, 'config.conf')
                with open(config_name, 'wt') as conf:
                    conf.write(f'scoring_function {sf}\n')
                    conf.write(f'search_speed {exhaustiveness}\n')
                    conf.write(f'protein_file {self.rec_path}\n')
                    conf.write(f'ligand_file {self.lig_path}\n')
                    conf.write(f'output_dir {join(out_dir, "results")}\n')
                    conf.write('write_multi_mol2 1\n')
                    conf.write(f'bindingsite_center {xc} {yc} {zc}\n')
                    conf.write(f'bindingsite_radius {r}\n')
                    conf.write(f'cluster_structures {self.num_poses}\n')
                    conf.write(f'cluster_rmsd {self.rmsd_tol}\n')
                config_paths.append(config_name)
        return config_paths

    def get_commands(self):
        config_paths = self.create_config()
        commands = []
        for config in config_paths:
            cmd = f'{self.exe} --mode screen {config}'
            commands.append(cmd)
        return commands, config_paths

    def run_commands(self):
        for i, command in enumerate(
                tqdm(self.commands, total=len(self.commands), unit='command')):
            path_split = self.config_paths[i].split(os.sep)
            log_base = os.sep.join(path_split[:-1])
            log_name = join(log_base, path_split[-2] + '.log')

            cmd_list = command.split()
            cmd_run = sp.Popen(self.bench + cmd_list, text=True,
                               stdout=sp.PIPE,
                               stderr=sp.PIPE)
            output, errors = cmd_run.communicate()
            cmn.write_string(output + errors, log_name)

    def get_info_dict(self):
        base_dict = cmn.recursive_defaultdict()

        rec_poses = cmn.recursive_finder('*bindingsite_fixed.mol2',
                                         self.out_dir)
        for rec_pose in rec_poses:
            path_split = rec_pose.split(os.sep)[:-2]
            id_ = path_split[-1]
            base_dict[id_]['rec_pose'] = rec_pose

        lig_poses = cmn.recursive_finder('docked_ligands.mol2', self.out_dir)
        for lig_pose in lig_poses:
            path_split = lig_pose.split(os.sep)[:-2]
            id_ = path_split[-1]
            base_dict[id_]['lig_pose'] = lig_pose

        scores = list(cmn.recursive_finder('ranking.csv', self.out_dir))
        for score in scores:
            path_split = score.split(os.sep)[:-2]
            id_ = path_split[-1]
            base_dict[id_]['score_pose'] = pd.read_csv(score).TOTAL_SCORE

        return base_dict

    def yield_filter_sort(self):
        info_dict = self.get_info_dict()
        for id_ in tqdm(info_dict, total=len(info_dict)):

            lig_parsed = cmn.Molecule(info_dict[id_]['lig_pose']).parse()
            rec_parsed = cmn.Molecule(info_dict[id_]['rec_pose']).parse()[0]
            rec_kdt = ckd(rec_parsed.getCoords())

            filtered_indices = []
            filtered_ligs = []
            for i, lig_ag in enumerate(lig_parsed):
                lig_kdt = ckd(lig_ag.hydrogen.getCoords())
                contacts = lig_kdt.query_ball_tree(rec_kdt, r=5)
                if any(contacts):
                    filtered_indices.append(i)
                    filtered_ligs.append(lig_ag)

            scores = info_dict[id_]['score_pose']
            out_dir = os.sep.join(
                info_dict[id_]['lig_pose'].split(os.sep)[:-2])
            scores.to_string(join(out_dir, 'scores.txt'))
            filtered_scores = scores[filtered_indices]
            filtered_scores.to_string(join(out_dir, 'scores_filtered.txt'))

            lig_ens = cmn.Molecule(info_dict[id_]['lig_pose']).get_ensemble()
            prd.writePDB(join(out_dir, 'poses.pdb'), lig_ens)

            lig_filtered = prd.Ensemble()
            [lig_filtered.addCoordset(x) for x in
             lig_ens.getCoordsets(filtered_indices)]
            lig_filtered.setAtoms(lig_ens[0].getAtoms())
            prd.writePDB(join(out_dir, 'poses_filtered.pdb'), lig_filtered)

# %%===========================================================================
# Debugging area
# =============================================================================
# from os.path import abspath
#
# root_dir = abspath('.')
# lig_mol2 = join(root_dir,
#                 'scripts/03_complexes_explorations/input_files/p97ND1/ligand.mol2')
# rec_mol2 = join(root_dir,
#                 'scripts/03_complexes_explorations/input_files/p97ND1/receptor.mol2')
#
# out_dir = join(root_dir, 'scripts/03_complexes_explorations/explorations/')
# n_poses = 2000
# rmsd_tol = 1.0
#
# plants_exe = '/home/roy.gonzalez-aleman/SoftWare/PLANTS1.2_64bit'
# plants_odir = join(out_dir, 'plants')
# plants_scores = ['plp', 'plp95', 'chemplp']
# plants_levels = ['speed1', 'speed2', 'speed4']
#
# self = Plants(plants_exe, rec_mol2, lig_mol2, n_poses, rmsd_tol,
#               plants_scores, plants_levels, plants_odir)
# self.run_commands()
# self.yield_filter_sort()
