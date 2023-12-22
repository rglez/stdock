# Created by roy.gonzalez-aleman at 23/11/2023
import os
import shutil
import subprocess as sp
from multiprocessing import Pool
from os.path import join, abspath, basename

import numpy as np
import prody as prd

import commons as cmn


def throw_error(cmd, label):
    if not cmd.returncode == 0:
        raise RuntimeError(f'Problems preparing {label.lower()} !')
    else:
        print(f'\n{label} has been correctly prepared')


class AutoDocker4:
    def __init__(self, lig_path, rec_path, ad4_path, adt_path, babel_path,
                 out_path, n_runs):

        # Parsed class arguments
        self.ga_run = 250
        self.ga_num_evals = 250000
        self.n_runs = n_runs
        self.rec_path = cmn.check_path(rec_path)
        self.lig_path = cmn.check_path(lig_path)
        self.ad4_path = cmn.check_path(ad4_path)
        self.adt_path = cmn.check_path(adt_path)
        self.out_path = abspath(cmn.check_path(out_path, check_exist=False))

        # Inferred from class arguments
        self.out_files = join(self.out_path, 'docking_files')
        self.out_poses = join(self.out_path, 'docking_poses')
        self.lig = prd.parsePDB(self.lig_path)
        self.rec = prd.parsePDB(self.rec_path)
        self.pythonsh = self.get_python_exe()
        self.autodock4 = join(self.ad4_path, 'autodock4')
        self.autogrid4 = join(self.ad4_path, 'autogrid4')
        self.babel_path = cmn.check_path(babel_path)

        # Common workflow
        self.build_hierarchy()
        self.prepare_ligand()
        self.prepare_receptor()
        # self.prepare_gpf()
        # self.compute_maps()
        # self.prepare_dpf()

    def get_python_exe(self):
        p = os.sep.join(self.adt_path.split(os.sep)[:-4] + ['bin', 'pythonsh'])
        return cmn.check_path(p)

    def build_hierarchy(self):
        os.makedirs(self.out_path)
        os.makedirs(self.out_files)
        os.makedirs(self.out_poses)

    def prepare_ligand(self):
        os.chdir(self.out_files)
        shutil.copy(self.lig_path, self.out_files)
        lig_script = next(cmn.recursive_finder('prep*lig*4.py', self.adt_path))
        cmd = sp.run([self.pythonsh, lig_script, '-l', self.lig_path])
        throw_error(cmd, 'Ligand')

    def prepare_receptor(self):
        os.chdir(self.out_files)
        shutil.copy(self.rec_path, self.out_files)
        rec_script = next(cmn.recursive_finder('pre*_rec*4.py', self.adt_path))
        cmd = sp.run([self.pythonsh, rec_script, '-r', self.rec_path])
        throw_error(cmd, 'Receptor')

    def define_dockbox(self):
        rec_coords = self.rec.getCoords()
        min_coords = rec_coords.min(axis=0)
        max_coords = rec_coords.max(axis=0)
        max_min = (max_coords - min_coords) + self.get_ligand_size()
        n_grids = (max_min / 0.375) // 2 * 2
        str_n_grids = ','.join([str(int(x)) for x in n_grids])
        rec_center = np.round(prd.calcCenter(self.rec), 1)
        center = ','.join([str(x) for x in rec_center])
        return center, str_n_grids

    def prepare_gpf(self):
        os.chdir(self.out_files)
        lig_qt = basename(self.lig_path) + 'qt'
        rec_qt = basename(self.rec_path) + 'qt'
        grid_center, n_pts = self.define_dockbox()
        gpf_script = next(cmn.recursive_finder('pre*_gpf*4.py', self.adt_path))
        cmd = sp.run(
            [self.pythonsh, gpf_script, '-l', lig_qt, '-r', rec_qt, '-p',
             f"npts={n_pts}", '-p', f"gridcenter={grid_center}"])
        throw_error(cmd, '.gpf')

    def get_ligand_size(self):
        coords = self.lig.getCoords()
        mini = coords.min(axis=0)
        maxi = coords.max(axis=0)
        return abs((maxi - mini).max())

    def compute_maps(self):
        os.chdir(self.out_files)
        rec_glg = basename(self.rec_path).split('.')[0] + '.glg'
        rec_gpf = basename(self.rec_path).split('.')[0] + '.gpf'
        print('\nComputing atomic affinity maps')
        cmd = sp.run([self.autogrid4, '-p', rec_gpf, '-l', rec_glg])
        throw_error(cmd, '.glg')

    def prepare_dpf(self):
        os.chdir(self.out_files)

        lig_qt = basename(self.lig_path) + 'qt'
        rec_qt = basename(self.rec_path) + 'qt'
        dpf_script = next(cmn.recursive_finder('pre*_dpf*4.py', self.adt_path))
        cmd = sp.run([self.pythonsh, dpf_script, '-l', lig_qt, '-r', rec_qt,
                      f'-p ga_num_evals={self.ga_num_evals}',
                      f'-p ga_run={self.ga_run}'],
                     stdout=sp.DEVNULL, stderr=sp.STDOUT)
        throw_error(cmd, '.dpf')

    def dock(self, i):
        os.chdir(self.out_files)
        dpf = next(cmn.recursive_finder('*dpf', '.'))
        cmd = sp.run([self.autodock4, '-p', dpf, '-l', f'run_{i}.dlg'],
                     stdout=sp.DEVNULL, stderr=sp.STDOUT)
        throw_error(cmd, '.dlg')

    def output_poses(self, run_id):
        os.chdir(self.out_files)
        dlg = next(cmn.recursive_finder(f'*_{run_id}.dlg', '.'))
        prefix = join(self.out_poses, f'poses_run_{run_id}')
        conf_script = next(
            cmn.recursive_finder('wr*_conf*dlg.py', self.adt_path))
        cmd = sp.run([self.pythonsh, conf_script, '-d', dlg, '-o', prefix],
                     stdout=sp.DEVNULL,
                     stderr=sp.STDOUT)
        if cmd.returncode != 0:
            raise RuntimeError('\nProblems extracting final poses')

    def output_scores(self):
        os.chdir(self.out_files)
        raw_dlg = list(cmn.recursive_finder('*dlg', '.'))
        dlg_files = sorted(raw_dlg, key=lambda x: int(
            basename(x).split('_')[1].split('.')[0]))

        flag = 'DOCKED: USER    Estimated Free Energy of Binding'
        scores = []
        for dlg in dlg_files:
            with open(dlg, 'rt') as out:
                run_scores = []
                for line in out:
                    if line.startswith(flag):
                        score = float(line.split('=')[1].split('kcal')[0])
                        run_scores.append(score)
            scores.extend(run_scores)

        scores_path = join(self.out_poses, 'poses_scores.txt')
        with open(scores_path, 'wt') as s:
            for score in scores:
                s.write(f'{score}\n')

    def convert_poses(self, run_id):
        os.chdir(self.out_poses)
        pdbqts = cmn.recursive_finder(f'*run_{run_id}*pdbqt', '.')
        print(f'\nConverting the run-{run_id} pdbqt files to pdb')
        for pdbqt in pdbqts:
            cmd = sp.run(
                [self.babel_path, '-ipdbqt', pdbqt, '-opdb', '-O', pdbqt[:-2]],
                stdout=sp.DEVNULL, stderr=sp.STDOUT)
            if cmd.returncode != 0:
                raise RuntimeError(
                    f'\nProblem converting {pdbqt} from pdbqt to pdb')
            os.remove(pdbqt)

    def create_pdb_traj(self):
        os.chdir(self.out_poses)
        pdb_files = list(cmn.recursive_finder('*.pdb', self.out_poses))
        sorted_pdbs = sorted(pdb_files,
                             key=lambda x:
                             (int(basename(x).split('_')[2]),
                              int(basename(x).split('_')[3].split('.')[0])))

        ensemble = prd.Ensemble()
        for i, pdb_ in enumerate(sorted_pdbs):
            ag_ = prd.parsePDB(pdb_)
            ensemble.addCoordset(ag_.getCoords())
            os.remove(pdb_)
        ensemble.setAtoms(ag_)
        prd.writePDB(join(self.out_poses, 'poses_coords.pdb'), ensemble)

    def reorder_traj(self):
        raw_traj_path = join(self.out_poses, 'poses_coords.pdb')
        scores_path = join(self.out_poses, 'poses_scores.txt')
        raw_traj = prd.parsePDB(raw_traj_path)
        raw_scores = pd.read_table(scores_path, header=None)[0].to_list()
        scores = np.asarray(raw_scores)
        order = scores.argsort()
        ordered_ensemble = prd.Ensemble()
        [ordered_ensemble.addCoordset(x) for x in
         raw_traj.getCoordsets()[order]]
        ordered_ensemble.setAtoms(raw_traj)
        prd.writePDB(join(self.out_poses, 'poses_coords_sorted.pdb'),
                     ordered_ensemble)

    def run(self):
        # Run parallel docking searches
        run_ids = list(range(1, self.n_runs + 1))
        pool = Pool()
        pool.map(self.dock, run_ids)

        # Convert poses into a pdb trajectory
        for run_id in run_ids:
            self.output_poses(run_id)
            self.convert_poses(run_id)
        self.create_pdb_traj()

        # Output scores
        self.output_scores()


# =============================================================================
#
# =============================================================================
root = '/home/roy.gonzalez-aleman/RoyHub/stdock/data/autodock/input_files'
ad4 = '/home/roy.gonzalez-aleman/SoftWare/autodock/x86_64Linux2/'
adt = '/home/roy.gonzalez-aleman/SoftWare/autodock/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/'

lig = '/home/roy.gonzalez-aleman/Posgrado_Bioinf/lig_Strytilcysteine.pdb'
rec = '/home/roy.gonzalez-aleman/Posgrado_Bioinf/rec_Eg5.pdb'
out_dir = '/home/roy.gonzalez-aleman/RoyHub/stdock/data/autodock/results'
babel = '/home/roy.gonzalez-aleman/miniconda3/bin/obabel'
self = AutoDocker4(lig, rec, ad4, adt, babel, out_dir, n_runs=8)


# self.run()

# %%
import pandas as pd
import matplotlib.pyplot as plt

scores_path = '/home/roy.gonzalez-aleman/RoyHub/stdock/data/autodock/results/docking_poses/poses_scores.txt'
scores = np.asarray(pd.read_table(scores_path, header=None)[0])
plt.hist(scores, bins=50, cumulative=False, alpha=0.75)
plt.show()

print(' '.join([str(x) for x in np.where(scores < -5)[0]]))

# %%
# =============================================================================
# Vina
# =============================================================================
from vina import Vina

receptor = '/home/roy.gonzalez-aleman/RoyHub/stdock/data/autodock/results/docking_files/5mcq.pdbqt'
ligand = '/home/roy.gonzalez-aleman/RoyHub/stdock/data/autodock/results/docking_files/ade.pdbqt'
center = [-44.2, 23.7, -33.1]
size = [78,  68, 60]

job = Vina()
job.set_receptor(receptor)
job.set_ligand_from_file(ligand)
job.compute_vina_maps(center=center, box_size=size)
job.dock(exhaustiveness=16, n_poses=5000)


job.energies(n_poses=5000)
job.write_poses('docking_results.pdbqt', overwrite=True, n_poses=10000, energy_range=1000)
