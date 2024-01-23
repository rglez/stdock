# Created by roy.gonzalez-aleman at 14/01/2024
"""
Prepares ligand and receptor molecules using Autodock Tools scripts
"""
import os
from collections import defaultdict
from os.path import basename, join, split

import numpy as np
import numpy_indexed as npi
import prody as prd

import commons as cmn
import proj_paths as pp


def get_ordering(lig_pdb_path, lig_pdbqt_path):
    """
    Gets the atomic ordering between a pdbqt file and its original pdb before
    conversion

    Args:
        lig_pdb_path: path to the ligand.pdb
        lig_pdbqt_path: path to the ligand.pdbqt

    Returns:
        Saves the ordering in the "qt_order.npy" file
    """
    parsed_pdb = prd.parsePDB(lig_pdb_path)
    parsed_qt = cmn.Molecule(lig_pdbqt_path).parse()[0]
    order = npi.indices(parsed_qt.getCoords(), parsed_pdb.getCoords())
    out_name = join(split(lig_pdb_path)[0], 'qt_order')
    np.save(out_name, order)


class Preparator:
    """
    Prototype of a Preparator class
    """

    def __init__(self, lig_path, rec_path, out_path, **kwargs):

        # Parse regular arguments
        self.lig_path = lig_path
        self.rec_path = rec_path
        self.out_path = out_path
        self.kwargs = kwargs

        # Prepare Ligand and Receptor
        self.prepare()

    def ligand_to_pdb(self):
        """
        Converts ligand to pdb if not already in this format

        Returns:
                path to the new ligand.pdb file if ligand in another format OR
                the same path to the original ligand.pdb
        """
        lig_base = basename(self.lig_path)
        lig_ext = lig_base.split('.')[-1]
        if lig_ext != 'pdb':
            lig_parsed = cmn.Molecule(self.lig_path).parse()[0]
            out_lig_pdb = join(self.out_path, f'{lig_base}.pdb')
            return prd.writePDB(out_lig_pdb, lig_parsed)
        else:
            return self.lig_path

    def check_kwarg_argument(self, arg_string):
        """
        Checks if the arg_string was passed to the kwargs defined on init

        Args:
            arg_string: argument string to be retrieved in the kwargs

        Returns:
            argument: the value of arg_string in kwargs. None if not defined
        """
        argument = self.kwargs.get(arg_string, None)
        if not argument:
            raise RuntimeError(f'The {arg_string} argument is mandatory')
        return argument

    def prepare_ligand(self):
        """
        Prepares the ligand
        """
        raise NotImplementedError

    def prepare_receptor(self):
        """
        Prepares the receptor
        """
        raise NotImplementedError

    def prepare(self):
        """
        Prepares ligand and receptor together with any other necessary step
        """
        raise NotImplementedError


class SporesPreparator(Preparator):
    """
    Prepares ligand and receptor molecules using SPORES scripts (PLANTS)
    """

    def prepare_mol(self, mol_path):
        """
        Prepares a molecule using SPORES

        Args:
            mol_path: a path to the molecule to be prepared
        """
        # Check spores
        spores_path = self.check_kwarg_argument('spores_path')
        cmn.check_path(spores_path)

        # Run spore
        mol_base_name = basename(mol_path).split('.')[0]
        mol_out_name = join(self.out_path, f'{mol_base_name}_spored.mol2')
        cmd = f'{spores_path} --mode settypes {mol_path} {mol_out_name}'
        cmn.shell_run(cmd.split())

    def prepare_receptor(self):
        self.prepare_mol(self.rec_path)
        pass

    def prepare_ligand(self):
        self.prepare_mol(self.lig_path)
        pass

    def prepare(self):
        self.prepare_receptor()
        self.prepare_ligand()
        pass


class AdtPreparator(Preparator):
    """
    Prepares ligand and receptor molecules using AutoDock Tools scripts
    """

    def prepare_ligand(self):
        # Check adt and pythonsh paths
        adt_path = self.check_kwarg_argument('adt_path')
        pythonsh = self.check_kwarg_argument('pythonsh')

        # Prepare ligand
        os.chdir(self.out_path)
        lig_prep_script = next(
            cmn.recursive_finder('prep*_lig*4.py', adt_path))
        lig_pdb_path = self.ligand_to_pdb()
        cmd_list = [pythonsh, lig_prep_script, '-l', lig_pdb_path, '-U', '\""']
        outputs, errors = cmn.shell_run(cmd_list)
        return outputs, errors

    def prepare_receptor(self):
        pass
        # Check adt and pythonsh paths
        adt_path = self.check_kwarg_argument('adt_path')
        pythonsh = self.check_kwarg_argument('pythonsh')

        # Prepare receptor
        os.chdir(self.out_path)
        rec_prep_script = next(
            cmn.recursive_finder('prep*_rec*4.py', adt_path))
        cmd_list = [pythonsh, rec_prep_script, '-r', self.rec_path]
        outputs, errors = cmn.shell_run(cmd_list)
        return outputs, errors

    def prepare(self):
        pass
        # Prepare ligand
        self.prepare_ligand()

        # Get order change from pdbqt conversion
        qt_name = basename(self.lig_path).split('.')[0] + '.pdbqt'
        lig_qt = join(self.out_path, qt_name)
        cmn.check_path(lig_qt)
        lig_pdb_path = self.ligand_to_pdb()
        get_ordering(lig_pdb_path, lig_qt)

        # Prepare receptor
        self.prepare_receptor()


# %%
# ==== Prepare folders hierarchy
proj_dir = pp.proj_dir
pythonsh = pp.pythonsh_exe
ad_tools_dir = pp.adtools_dir
spores_path = pp.spores_exe
all_files = join(proj_dir, 'data/external/coreset')

selected_log = join(proj_dir,
                    'scripts/01_complexes_selection/03_selection/report_selection.txt')
root_dir = split(selected_log)[0]
out_dir = join(proj_dir, 'scripts/02_complexes_preparation/01_prepared')
cmn.makedir_after_overwriting(out_dir)

# ==== Start the preparation
failed = defaultdict(list)

with open(selected_log, 'rt') as sele:
    for line in sele:
        # Get info from the selected cases log
        splitted = line.split(':')
        lig_path = splitted[1].strip()
        num_rots = int(splitted[0])

        # Create the preparation dir
        case = basename(lig_path).split('_')[1]
        case_dir_path = join(all_files, case)
        out_sub_dir = join(out_dir, case)
        cmn.makedir_after_overwriting(out_sub_dir)
        print(f'Preparing {case}')

        # Copy info into the preparation dir
        orig_pdb = join(case_dir_path, f'{case}_protein.pdb')
        rec_parsed = prd.parsePDB(orig_pdb).protein
        rec_pdb = join(out_sub_dir, f'{case}_protein.pdb')
        prd.writePDB(rec_pdb, rec_parsed)

        orig_pocket = join(case_dir_path, f'{case}_pocket.pdb')
        asr_resnums = np.unique(prd.parsePDB(orig_pocket).getResnums())
        np.save(join(out_sub_dir, 'bpocket'), asr_resnums)

        orig_lig = join(case_dir_path, f'{case}_ligand.mol2')
        lig_pdb = join(out_sub_dir, f'{case}_ligand.pdb')
        cmn.reformat_single_x2y(orig_lig, lig_pdb)

        # Prepare ligand & receptor with ADT
        try:
            adt = AdtPreparator(lig_pdb, rec_pdb, out_sub_dir,
                                adt_path=ad_tools_dir, pythonsh=pythonsh)
            print(f'{case} has been prepared with ADT')
        except RuntimeError:
            failed['adt'].append(case)
            print(f'\nPreparation of {case} with ADT failed')

        # Prepare ligand & receptor with SPORES
        try:
            self = SporesPreparator(lig_pdb, rec_pdb, out_sub_dir,
                                    spores_path=spores_path)
            print(f'{case} has been prepared with SPORE')
        except RuntimeError:
            failed['spore'].append(case)
            print(f'\nPreparation of {case} with SPORE failed')

        os.chdir(root_dir)
