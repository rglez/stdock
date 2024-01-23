# Created by roy.gonzalez-aleman at 03/01/2024
import os
import subprocess as sp
from os.path import join

from tqdm import tqdm

import commons as cmn
from programs.vina import Vina


class QvinaW(Vina):
    def get_commands(self):
        xc, yc, zc = self.get_rec_center()
        lig_size = cmn.get_longest_component(self.lig_parsed)
        xs, ys, zs = self.get_rec_axis() + lig_size
        commands = []
        for exhaustiveness in self.exhaustiveness_list:
            sub_dir = f'qvinaw_vina_{cmn.syn[exhaustiveness]}'
            out_dir = join(self.out_dir, sub_dir)
            cmd = (f'{self.exe}'
                   f' --receptor {self.rec_path}'
                   f' --ligand {self.lig_path}'
                   f' --center_x {xc}'
                   f' --center_y {yc}'
                   f' --center_z {zc}'
                   f' --size_x {xs}'
                   f' --size_y {ys}'
                   f' --size_z {zs}'
                   f' --seed 910623'
                   f' --exhaustiveness {exhaustiveness}'
                   f' --num_modes {self.num_poses}'
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
            cmn.write_string(output + errors, log_name)
