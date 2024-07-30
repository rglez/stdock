# Created by roy.gonzalez-aleman at 14/02/2024
import os
import shutil
import subprocess as sp
from os.path import join, split

import commons as cmn


def run_cmd(cmd):
    cmd_run = sp.Popen(cmd, text=False, stdout=sp.PIPE, stderr=sp.PIPE,
                       shell=True)
    return cmd_run.communicate()


class ADCP:
    def __init__(self, lig_pdb, rec_pdb, adfr_path, out_dir):
        # Create out dir & copy files
        self.out_dir = cmn.check_path(out_dir, check_exist=False)
        cmn.makedir_after_overwriting(self.out_dir)
        self.lig_pdb = shutil.copy(cmn.check_path(lig_pdb), self.out_dir)
        self.rec_pdb = shutil.copy(cmn.check_path(rec_pdb), self.out_dir)
        os.chdir(self.out_dir)

        # Check & parse ADFR executables
        self.adfr_path = cmn.check_path(adfr_path)
        self.reducer = cmn.check_path(join(adfr_path, 'reduce'))
        self.rec_prep = cmn.check_path(join(adfr_path, 'prepare_receptor'))
        self.lig_prep = cmn.check_path(join(adfr_path, 'prepare_ligand'))
        self.agfr = cmn.check_path(join(adfr_path, 'agfr'))


    def run(self):
        # 1. Reduce lig & rec
        lig_reduced = self.reduce(self.lig_pdb)
        rec_reduced = self.reduce(self.rec_pdb)

        # 2. Prepare lig & rec
        lig_qt = self.prepare_ligand(lig_reduced)
        rec_qt = self.prepare_receptor(rec_reduced)

        # 3. Generate maps
        self.generate_maps(lig_qt, rec_qt)

    def reduce(self, pdb_path):
        dirname, basename = split(pdb_path)
        out_pdb = join(dirname, f'{basename.split(".")[0]}_H.pdb')
        cmd = f'{self.reducer} {pdb_path} > {out_pdb}'
        stderr, stdout = run_cmd(cmd)
        return cmn.check_path(out_pdb)

    def prepare_ligand(self, lig_reduced):
        dirname, basename = split(lig_reduced)
        cmd = f'{self.lig_prep} -l {basename}'
        run_cmd(cmd)
        return cmn.check_path(f'{basename}qt')

    def prepare_receptor(self, rec_reduced):
        dirname, basename = split(rec_reduced)
        cmd = f'{self.rec_prep} -r {basename}'
        run_cmd(cmd)
        return cmn.check_path(f'{basename}qt')

    def generate_maps(self, lig_qt, rec_qt):
        cmd = f'{self.agfr} -r {rec_qt} -l {lig_qt} -o complex'
        return run_cmd(cmd)


# agfr -r 5GRD_recH.pdbqt -l 5GRD_pepH.pdbqt -o 5GRD
# =============================================================================
# Debugging area
# =============================================================================
# adfr_path = '/home/roy.gonzalez-aleman/ADFRsuite-1.0/bin/'
# rec_pdb = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/STDCK2a-300v2/ck2_alpha.pdb'
# lig_pdb = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/STDCK2a-300v2/CIGB300_MINIM.pdb'
# out_dir = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/ADCP'
# self = ADCP(lig_pdb, rec_pdb, adfr_path, out_dir)
# self.run()
