# Created by roy.gonzalez-aleman at 14/01/2024
import os
import shutil
import subprocess as sp
from os.path import basename, join, split, abspath

import numpy as np
import numpy_indexed as npi
import prody as prd

import commons as cmn
import root
# todo: be careful with paths management here

def throw_error(cmd, label):
    if not cmd.returncode == 0:
        raise RuntimeError(f'Problems preparing {label.lower()} !')
    else:
        print(f'\n{label} has been correctly prepared')


def get_ordering(lig_pdb_path, lig_pdbqt_path):
    parsed_pdb = prd.parsePDB(lig_pdb_path)
    parsed_qt = root.Molecule(lig_pdbqt_path).parse()[0]
    order = npi.indices(parsed_qt.getCoords(), parsed_pdb.getCoords())
    out_name = join(split(lig_pdb_path)[0], 'qt_order')
    np.save(out_name, order)


class Preparator:
    def __init__(self, pythonsh, adt_path, rec_path, lig_path, out_path):
        # Parse arguments
        self.pythonsh = cmn.check_path(pythonsh)
        self.adt_path = cmn.check_path(adt_path)
        self.rec_path = cmn.check_path(rec_path)
        self.lig_path = cmn.check_path(lig_path)
        self.lig_pdb_path = self.ligand_to_pdb()
        self.out_path = cmn.check_path(out_path)

        # Find adt preparation scripts
        self.lig_prep_script = next(
            cmn.recursive_finder('prep*_lig*4.py', adt_path))
        self.rec_prep_script = next(
            cmn.recursive_finder('prep*_rec*4.py', adt_path))

        # Trigger preparations
        self.prepare()

    def ligand_to_pdb(self):
        lig_ext = basename(self.lig_path).split('.')[-1]
        if lig_ext != 'pdb':
            lig_parsed = root.Molecule(self.lig_path).parse()[0]
            return prd.writePDB(f'{self.lig_path}.pdb', lig_parsed)
        else:
            return self.lig_path

    def prepare_ligand(self):
        # Prepare ligand
        os.chdir(self.out_path)
        cmd = sp.run(
            [self.pythonsh, self.lig_prep_script, '-l', self.lig_path, '-U',
             '\""'])
        throw_error(cmd, 'Ligand')

    def prepare_receptor(self):
        os.chdir(self.out_path)
        cmd = sp.run([self.pythonsh, self.rec_prep_script, '-r', rec_path])
        throw_error(cmd, 'Receptor')

    def prepare(self):
        # Prepare ligand
        self.prepare_ligand()

        # Get order change from pdbqt conversion
        qt_name = basename(self.lig_path).split('.')[0] + '.pdbqt'
        lig_qt = join(self.out_path, qt_name)
        cmn.check_path(lig_qt)
        get_ordering(self.lig_pdb_path, lig_qt)

        # Prepare receptor
        self.prepare_receptor()


# =============================================================================
# User-defined parameters
# =============================================================================
pythonsh = '/home/roy.gonzalez-aleman/SoftWare/autodock/mgltools_x86_64Linux2_1.5.7/bin/pythonsh'
ad_tools_dir = '/home/roy.gonzalez-aleman/SoftWare/autodock/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/'

all_files = '/home/roy.gonzalez-aleman/DataHub/data_raw/v2013-core/'
selected_log = abspath('scripts/01_complexes_selection/03_selection/report_selection.txt')
root_dir = split(abspath(selected_log))[0]
out_dir = abspath('scripts/02_complexes_preparation/01_prepared_with_adt')
failed = []
with open(selected_log, 'rt') as sele:
    for line in sele:
        # Get info from the selected cases log
        splitted = line.split(':')
        lig_path = splitted[1].strip()
        num_rots = int(splitted[0])

        # Get all necessary file paths
        case = basename(lig_path).split('_')[1]
        case_dir_path = join(all_files, case)
        rec_path = join(case_dir_path, f'{case}_protein.pdb')
        rec_parsed = prd.parsePDB(rec_path).protein
        prd.writePDB(rec_path, rec_parsed)
        lig_path = join(case_dir_path, f'{case}_ligand.mol2')
        pocket_path = join(case_dir_path, f'{case}_pocket.pdb')

        # Copy info to the output folder
        out_sub_dir = join(out_dir, case)
        os.makedirs(out_sub_dir, exist_ok=True)
        shutil.copy(rec_path, out_sub_dir)
        shutil.copy(lig_path, out_sub_dir)
        asr_resnums = np.unique(prd.parsePDB(pocket_path).getResnums())
        np.save(join(out_sub_dir, 'bpocket'), asr_resnums)

        # Prepare ligand & receptor
        try:
            self = Preparator(pythonsh, ad_tools_dir, rec_path, lig_path,
                              out_sub_dir)
        except RuntimeError:
            failed.append(case)
        os.chdir(root_dir)
