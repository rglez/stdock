# Created by roy.gonzalez-aleman at 02/01/2024
import os
import re
import subprocess as sp
from os.path import join, split

import prody as prd
from scipy.spatial import cKDTree as ckd
from tqdm import tqdm

import commons as cmn
import root


class Vina(root.Program):
    def get_commands(self):
        xc, yc, zc = self.get_rec_center()
        lig_size = root.get_longest_component(self.lig_parsed)
        xs, ys, zs = self.get_rec_axis() + lig_size
        commands = []
        for exhaustiveness in self.exhaustiveness_list:
            for sf in self.scoring_functions:
                sub_dir = f'vina_{sf}_{root.syn[exhaustiveness]}'
                out_dir = join(self.out_dir, sub_dir)
                cmd = (f'{self.exe}'
                       f' --receptor {self.rec_path}'
                       f' --ligand {self.lig_path}'
                       f' --scoring {sf}'
                       f' --center_x {xc}'
                       f' --center_y {yc}'
                       f' --center_z {zc}'
                       f' --size_x {xs}'
                       f' --size_y {ys}'
                       f' --size_z {zs}'
                       f' --seed 910623'
                       f' --exhaustiveness {exhaustiveness}'
                       f' --num_modes {self.num_poses}'
                       f' --min_rmsd {self.rmsd_tol}'
                       f' --energy_range 10000'
                       f' --out {join(out_dir, sub_dir)}.pdbqt')
                commands.append(cmd)
        return commands, None

    def run_commands(self):
        for i, command in enumerate(
                tqdm(self.commands, total=len(self.commands), unit='command')):
            odir = os.sep.join(
                command.split('--out')[-1].strip().split(os.sep)[:-1])
            os.makedirs(odir, exist_ok=True)

            log_name = join(odir, odir.split(os.sep)[-1] + '.log')

            cmd_list = command.split()
            cmd_run = sp.Popen(self.bench + cmd_list, text=True,
                               stdout=sp.PIPE,
                               stderr=sp.PIPE)
            output, errors = cmd_run.communicate()
            root.write_string(output + errors, log_name)

    def yield_filter_sort(self):
        string_score = 'REMARK VINA RESULT:.*'

        for pdbqt in cmn.recursive_finder('*pdbqt', self.out_dir):
            outdir = split(pdbqt)[0]

            # Get all scores
            with open(pdbqt, 'rt') as inp:
                file_string = inp.read()
            score_raw = re.findall(string_score, file_string)
            scores = [x.split(':')[1].split()[0] for x in score_raw]
            scores_path = join(outdir, 'scores.txt')
            with open(scores_path, 'wt') as out:
                for i, x in enumerate(scores):
                    out.write(f'{i}    {x}\n')

            # Get all poses
            ensemble = root.Molecule(pdbqt).get_ensemble()
            out_name = join(outdir, 'poses.pdb')
            prd.writePDB(out_name, ensemble)

            # Get filtered indices
            rec_kdt = ckd(self.rec_parsed.getCoords())
            lig_parsed = root.Molecule(pdbqt).parse()
            filtered_indices, filtered_ligs = root.get_filtered_indices(
                rec_kdt, lig_parsed)

            # Get filtered poses
            lig_filtered = prd.Ensemble()
            [lig_filtered.addCoordset(x.getCoords()) for x in filtered_ligs]
            lig_filtered.setAtoms(filtered_ligs[0])
            prd.writePDB(join(outdir, 'poses_filtered.pdb'), lig_filtered)

            # Get filtered scores
            with open(join(outdir, 'scores_filtered.txt'), 'wt') as out:
                for index in filtered_indices:
                    out.write(f'{index}    {scores[index]}\n')

# =============================================================================
# Debugging area
# =============================================================================

# inputs_dir = 'scripts/03_complexes_explorations/input_files/'
# cases = os.listdir(inputs_dir)
# for case in cases:
#     lig_qt = join(inputs_dir, case, 'ligand.pdbqt')
#     rec_qt = join(inputs_dir, case, 'receptor.pdbqt')
#     out_dir = f'scripts/03_complexes_explorations/explorations/{case}'
#
# python_exe = '/home/roy.gonzalez-aleman/miniconda3/envs/stdock/bin/python'
# n_poses = 2000
# rmsd_tol = 1.0
#
# vina_exe = "/home/roy.gonzalez-aleman/SoftWare/vina_1.2.5_linux_x86_64"
# vina_odir = join(out_dir, 'vina')
# vina_scores = ['ad4', 'vina', 'vinardo']
# vina_levels = [800, 80, 8]
# self = Vina(vina_exe, rec_qt, lig_qt, n_poses, rmsd_tol,
#             vina_scores, vina_levels, vina_odir)
