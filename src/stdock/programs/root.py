# Created by roy.gonzalez-aleman at 01/01/2024
import os
import tempfile
from os.path import split, join

import numpy as np
import prody as prd
from openbabel import pybel
from scipy.spatial import cKDTree as ckd

import commons as cmn

syn = {8: 'low',
       80: 'medium',
       800: 'high',
       'speed1': 'high',
       'speed2': 'medium',
       'speed4': 'low',
       250: 'low',
       2500: 'medium',
       25000: 'high'}


def get_filtered_indices(rec_kdt, lig_parsed):
    filtered_indices = []
    filtered_ligs = []
    for i, lig_ag in enumerate(lig_parsed):
        try:
            lig_kdt = ckd(lig_ag.getCoords())
            contacts = lig_kdt.query_ball_tree(rec_kdt, r=5)
            if any(contacts):
                filtered_indices.append(i)
                filtered_ligs.append(lig_ag)
        except AttributeError:
            print(f'WARNING: Frame {i} not filtered')
    return filtered_indices, filtered_ligs


def write_string(string, path):
    with open(path, 'wt') as out:
        out.write(string)
    return path


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
        rec_obj = pybel.readfile(extension, self.mol_path)

        with tempfile.TemporaryDirectory() as d:
            parsed = []
            for i, obj in enumerate(rec_obj):
                temp_mol = join(d, f'temp_molecule_{i}.pdb')
                obj.write("pdb", temp_mol)
                parsed.append(prd.parsePDB(temp_mol))
        return parsed

    def get_ensemble(self):
        parsed = self.parse()
        ensemble = prd.Ensemble()
        for i, ag in enumerate(parsed):
            try:
                ensemble.addCoordset(ag.getCoords())
            except AttributeError:
                print(f'WARNING: Frame {i} not parsed.')
                pass
        ensemble.setAtoms(parsed[0])
        return ensemble


class Program:
    def __init__(self, exe_path, rec_path, lig_path, n_poses, rmsd_tol,
                 scoring_functions, exhaustiveness_list, odir):
        # Program executable path
        self.exe = cmn.check_path(exe_path)

        # Receptor & Ligand parsing
        self.rec_path = cmn.check_path(rec_path)
        self.rec_parsed = Molecule(self.rec_path).parse()[0]
        self.lig_path = cmn.check_path(lig_path)
        self.lig_parsed = Molecule(self.lig_path).parse()[0]

        # Gather all scoring functions
        self.scoring_functions = scoring_functions

        # Exhaustiveness
        self.exhaustiveness_list = exhaustiveness_list

        # Check numeric arguments
        self.num_poses = cmn.check_numeric_in_range('n_poses', n_poses, int, 1,
                                                    cmn.inf_int)
        self.rmsd_tol = cmn.check_numeric_in_range('rmsd_tol', rmsd_tol, float,
                                                   0.1, cmn.inf_float)
        # Create output dir
        self.out_dir = odir
        os.makedirs(self.out_dir, exist_ok=True)

        # Standardize programs' running
        self.commands, self.config_paths = self.get_commands()
        self.bench = ['/usr/bin/time', '-v']

    def get_rec_axis(self, output_min_max=False):
        rec_coords = self.rec_parsed.getCoords()
        min_coords = rec_coords.min(axis=0)
        max_coords = rec_coords.max(axis=0)
        if output_min_max:
            return max_coords - min_coords, min_coords, max_coords
        return max_coords - min_coords

    def get_rec_center(self):
        rec_parsed = Molecule(self.rec_path).parse()[0]
        return np.round(prd.calcCenter(rec_parsed), 1)

    def get_commands(self):
        raise NotImplementedError

    def run_commands(self):
        raise NotImplementedError

    def yield_filter_sort(self):
        raise NotImplementedError
