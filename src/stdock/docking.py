# Created by roy.gonzalez-aleman at 26/12/2023
import os
import tempfile
from os.path import split, join

import numpy as np
import prody as prd
from openbabel import pybel

import commons as cmn

syn = {8: 'low',
       80: 'medium',
       800: 'high',
       'speed1': 'high',
       'speed2': 'medium',
       'speed4': 'low'}


def get_longest_component(parsed_mol):
    coords = parsed_mol.getCoords()
    mini = coords.min(axis=0)
    maxi = coords.max(axis=0)
    return abs((maxi - mini).max())


class Molecule:
    def __init__(self, molecule_path):
        self.mol_path = molecule_path
        cmn.check_path(self.mol_path)

    def parse(self):
        root, base = split(self.mol_path)
        extension = base.split('.')[-1]
        rec_obj = next(pybel.readfile(extension, self.mol_path))

        with tempfile.TemporaryDirectory() as d:
            temp_mol = join(d, 'temp_molecule.pdb')
            rec_obj.write("pdb", temp_mol)
            parsed = prd.parsePDB(temp_mol)
        return parsed


class Program:
    def __init__(self, exe_path, rec_path, lig_path, n_poses, rmsd_tol,
                 scoring_functions, exhaustiveness_list, odir):
        # Program executable path
        self.exe = cmn.check_path(exe_path)

        # Receptor & Ligand parsing
        self.rec_path = cmn.check_path(rec_path)
        self.rec_parsed = Molecule(self.rec_path).parse()

        self.lig_path = cmn.check_path(lig_path)
        self.lig_parsed = Molecule(self.lig_path).parse()

        # Gather all scoring functions
        self.scoring_functions = scoring_functions

        # Exhaustiveness
        self.exhaustiveness_list = exhaustiveness_list

        # Check numeric arguments & create output dir
        self.num_poses = cmn.check_numeric_in_range('n_poses', n_poses, int, 1,
                                                    cmn.inf_int)
        self.rmsd_tol = cmn.check_numeric_in_range('rmsd_tol', rmsd_tol, float,
                                                   0.1, cmn.inf_float)
        self.out_dir = odir
        os.makedirs(self.out_dir, exist_ok=True)

    def get_rec_axis(self):
        rec_coords = self.rec_parsed.getCoords()
        min_coords = rec_coords.min(axis=0)
        max_coords = rec_coords.max(axis=0)
        return max_coords - min_coords

    def get_rec_center(self):
        rec_parsed = Molecule(self.rec_path).parse()
        return np.round(prd.calcCenter(rec_parsed), 1)

    def get_commands(self):
        raise NotImplementedError


class QvinaW(Program):
    def get_commands(self):
        xc, yc, zc = self.get_rec_center()
        lig_size = get_longest_component(self.lig_parsed)
        xs, ys, zs = self.get_rec_axis() + lig_size
        commands = []
        for exhaustiveness in self.exhaustiveness_list:
            out_dir = join(self.out_dir,
                           f'qvinaw_vina_{syn[exhaustiveness]}')
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
                   f' --out {out_dir}')
            commands.append(cmd)
        return commands


class Vina(Program):
    def get_commands(self):
        xc, yc, zc = self.get_rec_center()
        lig_size = get_longest_component(self.lig_parsed)
        xs, ys, zs = self.get_rec_axis() + lig_size
        commands = []
        for exhaustiveness in self.exhaustiveness_list:
            for sf in self.scoring_functions:
                out_dir = join(self.out_dir,
                               f'vina_{sf}_{syn[exhaustiveness]}')
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
                       f' --out {out_dir}')
                commands.append(cmd)
        return commands


class Gnina(Program):
    def get_commands(self):
        commands = []
        for exhaustiveness in self.exhaustiveness_list:
            for sf in self.scoring_functions:
                out_dir = join(self.out_dir,
                               f'gnina_{sf}_{syn[exhaustiveness]}')
                cmd = (f'{self.exe}'
                       f' --receptor {self.rec_path}'
                       f' --ligand {self.lig_path}'
                       f' --scoring {sf}'
                       f' --autobox_ligand {self.rec_path}'
                       f' --cnn_scoring rescore'
                       f' --seed 910623'
                       f' --exhaustiveness {exhaustiveness}'
                       f' --num_modes {self.num_poses}'
                       f' --min_rmsd_filter {self.rmsd_tol}'
                       f' --out {out_dir}'
                       )
                commands.append(cmd)
        return commands


class Smina(Program):
    def get_commands(self):
        commands = []
        for exhaustiveness in self.exhaustiveness_list:
            for sf in self.scoring_functions:
                out_dir = join(self.out_dir,
                               f'smina_{sf}_{syn[exhaustiveness]}')
                cmd = (f'{self.exe}'
                       f' --receptor {self.rec_path}'
                       f' --ligand {self.lig_path}'
                       f' --scoring {sf}'
                       f' --autobox_ligand {self.rec_path}'
                       f' --seed 910623'
                       f' --exhaustiveness {exhaustiveness}'
                       f' --num_modes {self.num_poses}'
                       f' --min_rmsd_filter {self.rmsd_tol}'
                       f' --energy_range 10000'
                       f' --out {out_dir}'
                       f' --log {out_dir}'
                       )
                commands.append(cmd)
        return commands


class Plants(Program):
    def create_config(self):
        xc, yc, zc = self.get_rec_center()
        lig_size = get_longest_component(self.lig_parsed)
        r = self.get_rec_axis().max() / 2 + lig_size

        config_paths = []
        for exhaustiveness in self.exhaustiveness_list:
            for sf in self.scoring_functions:
                suffix = f'plants_{sf}_{syn[exhaustiveness]}'
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
        commands = []
        for config in self.create_config():
            cmd = f'{self.exe} --mode screen {config}'
            commands.append(cmd)
        return commands
