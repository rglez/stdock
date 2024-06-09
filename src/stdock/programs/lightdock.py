# Created by roy.gonzalez-aleman at 06/02/2024

import multiprocessing
import os
import shutil
import subprocess as sp
from os.path import join, basename, dirname

import numpy as np
import prody as prd
from prody import parsePDB

import commons as cmn


def run_command(command_string, log_name=None):
    """
    Run a command string using subprocess

    Args:
        command_string: a string of a valid command to execute
        log_name: a path to a log file for saving command's outputs
    """
    cmd_run = sp.Popen(command_string.split(), text=True, stdout=sp.PIPE,
                       stderr=sp.PIPE)
    output, errors = cmd_run.communicate()
    if log_name:
        cmn.write_string(output + errors, log_name)


def parse_pdb(pdb_path):
    structure = parsePDB(pdb_path, chain='L')
    structure.setTitle(pdb_path)
    yield structure


def swarm_to_dcd(swarm_dir):
    """
    Convert all PDB files inside a swarm dir into a DCD file

    Args:
        swarm_dir: path to a swarm dir containing PDB files
    """

    pdbs_raw = cmn.recursive_finder('*pdb', swarm_dir)
    pdb_function = lambda x: int(basename(x).split('.')[0].split('_')[-1])
    pdb_sorted = sorted(pdbs_raw, key=pdb_function)

    ensemble = prd.Ensemble()
    for pdb_path in pdb_sorted:
        pdb_parsed = prd.parsePDB(pdb_path, chain='L')
        pdb_coords = pdb_parsed.getCoords()
        ensemble.addCoordset(pdb_coords)
    ensemble.setCoords(pdb_coords)

    pdb_dirname = dirname(pdb_path)
    out_dcd = f'{basename(pdb_dirname)}.dcd'
    prd.writeDCD(join(pdb_dirname, out_dcd), ensemble)
    [os.remove(x) for x in pdb_sorted]


class LightDock:
    """
    LightDock software manager
    """

    def __init__(self, rec_pdb, lig_pdb, out_dir, num_steps=1,
                 num_swarms=100, num_glowworms=5, sf='tobi'):
        """
        Args:
            rec_pdb: path to the receptor pdb file
            lig_pdb: path to the ligand pdb file
            out_dir: path to the output dir
            num_steps: number of docking steps per swarm
            sf: scoring function to use
        """

        # Parse arguments
        self.num_steps = num_steps
        self.num_swarms = num_swarms
        self.num_glowworms = num_glowworms
        self.sf = sf
        self.out_dir = cmn.makedir_after_overwriting(out_dir)
        self.out_dir = out_dir
        os.chdir(self.out_dir)

        # Check rec, lig and rename their chains to R and L
        self.rec_path = self.check_and_chain_rename(rec_pdb, 'R')
        self.lig_path = self.check_and_chain_rename(lig_pdb, 'L')

        # Public attributes to compute by run method
        self._out_files = None
        self.light_dock_scores = None
        self.poses_dcd = None

        # Run workflow
        self.run()

    def check_and_chain_rename(self, pdb_path, chain):
        """
        Check existence of passed pdb files and assign a chain ID
        Args:
            pdb_path: path to a pdb file
            chain: chain ID to assign to the passed pdb file

        Returns:
            the path to the pdb file with a new chain ID set
        """
        cmn.check_path(pdb_path)
        pdb_path_raw = shutil.copy(pdb_path, self.out_dir)
        pdb_renamed = prd.parsePDB(pdb_path_raw)
        pdb_renamed.setChids(chain)
        return prd.writePDB(pdb_path_raw, pdb_renamed)

    def run(self):
        """
        Run LightDock workflow
        """
        # 1. Setup
        print("Defining LightDock's setup")
        self.setup()

        # 2. Dock
        print(f'Starting LightDock: {self.num_swarms} swarm(s),'
              f' {self.num_glowworms} glowworm(s), and'
              f' {self.num_steps} step(s) using {self.sf}')
        self.dock()

        # 3. Score
        print("Parsing LightDock's scores from out files")
        self.score()

        # 4. Extract
        print('Extracting LightDock poses as PDB files')
        self.extract()

        # 5. Convert PDB poses to DCD trajectory
        print('Converting LightDock poses to DCD files')
        self.swarm_conversion()
        self.poses_dcd = self.concat_swarms()

    def setup(self):
        """
        Set up the job files
        """
        setup = shutil.which('lightdock3_setup.py')
        log_name = join(self.out_dir, 'setup.log')
        # cmd = (f'{setup} {self.rec_path} {self.lig_path} --noxt --now -r /home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/cigb300_restraints.txt'
        #        f' -g {self.num_glowworms} -s {self.num_swarms}')
        cmd = (f'{setup} {self.rec_path} {self.lig_path} --noxt --now -r /home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/cigb300_restraints.txt'
               f' -g {self.num_glowworms} -s {self.num_swarms}')
        run_command(cmd, log_name)

    def dock(self):
        """
        Dock
        """
        docker = shutil.which('lightdock3.py')
        log_name = join(self.out_dir, 'dock.log')
        cmd = f'{docker} setup.json {self.num_steps} -s {self.sf} -min'
        run_command(cmd, log_name)

    def score(self):
        """
        Parse poses score from the last gso*.out file produced by LightDock
        """
        # Get out files and export for future usage
        out_raw = list(cmn.recursive_finder('gso_*.out', self.out_dir))
        _ = [int(basename(x).split('.')[0].split('_')[-1]) for x in out_raw]
        num = max(_)
        out_files = list(cmn.recursive_finder(f'gso_{num}.out', self.out_dir))
        self._out_files = sorted(out_files, key=lambda x: int(
            x.split(os.sep)[-2].split('_')[-1].split('.')[0]))

        # Parse scores from out files
        scores = []
        for out_file in self._out_files:
            with open(out_file, 'rt') as of:
                of.readline()
                for line in of:
                    scores.append(float(line.split()[-1]))
        self.light_dock_scores = np.asarray(scores)

    def extract(self):
        """
        Extract all generated poses as PDB files
        """
        extractor = shutil.which('lgd_generate_conformations.py')
        out_files = self._out_files

        extract_commands = []
        for gso in out_files:
            cmd_extract = f'{extractor} {self.rec_path} {self.lig_path} {gso} {cmn.inf_int}'
            extract_commands.append(cmd_extract)

        pool = multiprocessing.Pool(os.cpu_count())
        pool.map(run_command, extract_commands)

    def get_swarm_dirs(self):
        """
        Get swarms dirs path

        Returns:
            list of absolute paths to swarms dirs
        """
        swarm_dirs_relative = filter(lambda x: x.startswith('swarm'),
                                     os.listdir(self.out_dir))
        return [join(self.out_dir, x) for x in swarm_dirs_relative]

    def swarm_conversion(self):
        """
        Parallel conversion of all PDB files inside swarm dirs into DCD files
        """
        swarm_dirs = self.get_swarm_dirs()
        pool = multiprocessing.Pool(os.cpu_count())
        pool.map(swarm_to_dcd, swarm_dirs)
        pool.close()
        pool.join()

    def concat_swarms(self):
        """
        Concatenate all DCD files inside swarm dirs into a single DCD

        Returns:
            path to the concatenated DCD file
        """
        swarms_dcd = cmn.recursive_finder('*.dcd', self.out_dir)
        swarms_sorted = sorted(swarms_dcd, key=lambda x: int(
            basename(x).split('.')[0].split('_')[-1]))
        traj = prd.Trajectory(swarms_sorted[0])
        [traj.addFile(x) for x in swarms_sorted[1:]]
        out_name = join(self.out_dir, 'swarms')
        swarm_dcd_path = prd.writeDCD(f'{out_name}.dcd', traj)
        [os.remove(x) for x in swarms_dcd]
        return swarm_dcd_path


# =============================================================================
# Debugging area
# =============================================================================
# lig_pdb = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/STDCK2a-300v2/CIGB300_MINIM.pdb'
# rec_pdb = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/STDCK2a-300v2/ck2_alpha.pdb'
# out_dir = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/stdock-CIGB300/DOCKING'
# sf = 'tobi'
# sf = 'pisa'
# sf = 'mj3h'
# num_steps = 1
# num_swarms = 8
# num_glowworms = 100
#
# self = LightDock(rec_pdb, lig_pdb, out_dir, num_steps=num_steps,
#                  num_glowworms=num_glowworms, num_swarms=num_swarms, sf=sf)
