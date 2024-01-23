# Created by roy.gonzalez-aleman at 03/01/2024
from os.path import join

import commons as cmn
from programs.smina import Smina


class Gnina(Smina):
    def get_commands(self):
        commands = []
        for exhaustiveness in self.exhaustiveness_list:
            for sf in self.scoring_functions:
                sub_dir = f'gnina_{sf}_{cmn.syn[exhaustiveness]}'
                out_dir = join(self.out_dir, sub_dir)
                cmd = (f'{self.exe}'
                       f' --receptor {self.rec_path}'
                       f' --ligand {self.lig_path}'
                       f' --scoring {sf}'
                       f' --autobox_ligand {self.rec_path}'
                       f' --cnn_scoring none'
                       f' --seed 910623'
                       f' --exhaustiveness {exhaustiveness}'
                       f' --num_modes {self.num_poses}'
                       f' --min_rmsd_filter {self.rmsd_tol}'
                       f' --out {join(out_dir, sub_dir)}.pdbqt'
                       )
                commands.append(cmd)
        return commands, None
