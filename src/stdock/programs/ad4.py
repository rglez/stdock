# Created by roy.gonzalez-aleman at 27/12/2023
import argparse
import os
import shutil
import subprocess as sp
from os.path import basename, join, split

import numpy as np
import numpy_indexed as npi
import prody as prd
from scipy.spatial import cKDTree as ckd

import commons as cmn
from programs import root
from programs.root import Program, syn


def throw_error(cmd, label):
    if not cmd.returncode == 0:
        raise RuntimeError(f'Problems preparing {label.lower()} !')
    else:
        print(f'\n{label} has been correctly prepared')


def get_ordering(lig_pdb_path):
    lig_pdb_qt_path = f'{lig_pdb_path}qt'
    parsed_pdb = prd.parsePDB(lig_pdb_path)
    parsed_qt = root.Molecule(lig_pdb_qt_path).parse()[0]
    order = npi.indices(parsed_qt.getCoords(), parsed_pdb.getCoords())
    out_name = join(split(lig_pdb_path)[0], 'qt_order')
    np.save(out_name, order)


class AutoDock4(Program):
    def __init__(self, ad4_path, adt_path, babel_path, exe_path, rec_path,
                 lig_path, n_poses, rmsd_tol, scoring_functions,
                 exhaustiveness_list, odir):
        super().__init__(exe_path, rec_path, lig_path, n_poses, rmsd_tol,
                         scoring_functions, exhaustiveness_list, odir)

        # Parsed class arguments
        self.n_poses = n_poses
        self.exhaustiveness_list = exhaustiveness_list

        self.rec_path = cmn.check_path(rec_path)
        self.lig_path = cmn.check_path(lig_path)
        self.ad4_path = cmn.check_path(ad4_path)
        self.adt_path = cmn.check_path(adt_path)
        self.out_path = join(odir, 'autodock4',
                             f'autodock4_{syn[self.exhaustiveness_list]}')

        # # Inferred from class arguments
        self.out_files = join(self.out_path, 'docking_files')
        self.out_poses = join(self.out_path, 'docking_poses')
        self.lig = prd.parsePDB(self.lig_path)
        self.rec = prd.parsePDB(self.rec_path)
        self.pythonsh = self.get_python_exe()
        self.autodock4 = join(self.ad4_path, 'autodock4')
        self.autogrid4 = join(self.ad4_path, 'autogrid4')
        self.babel_path = cmn.check_path(babel_path)

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
        get_ordering(self.lig_path)

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
                      f'-p ga_num_evals={self.exhaustiveness_list}',
                      f'-p ga_run={self.n_poses}'],
                     stdout=sp.DEVNULL, stderr=sp.STDOUT)
        throw_error(cmd, '.dpf')

    def dock(self):
        os.chdir(self.out_files)
        dpf = next(cmn.recursive_finder('*dpf', '.'))
        cmd = sp.run([self.autodock4, '-p', dpf, '-l', f'run.dlg'],
                     stdout=sp.DEVNULL, stderr=sp.STDOUT)
        throw_error(cmd, '.dlg')

    def output_poses(self):
        os.chdir(self.out_files)
        dlg = next(cmn.recursive_finder(f'run.dlg', '.'))
        prefix = join(self.out_poses, f'poses_run')
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
        dlg_files = sorted(raw_dlg)

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

    def convert_poses(self):
        os.chdir(self.out_poses)
        pdbqts = cmn.recursive_finder(f'*run*pdbqt', '.')
        print(f'\nConverting the pdbqt files to pdb')
        for pdbqt in pdbqts:
            cmd = sp.run(
                [self.babel_path, '-ipdbqt', pdbqt, '-opdb', '-O', pdbqt[:-2]],
                stdout=sp.DEVNULL, stderr=sp.STDOUT)
            if cmd.returncode != 0:
                raise RuntimeError(
                    f'\nProblem converting {pdbqt} from pdbqt to pdb')
            os.remove(pdbqt)

    def create_pdbqt(self):
        os.chdir(self.out_poses)
        pdb_files = list(cmn.recursive_finder('*.pdb', self.out_poses))
        sorted_pdbs = sorted(pdb_files,
                             key=lambda x: int(
                                 basename(x).split('_')[2].split('.')[0]))

        ensemble = prd.Ensemble()
        for i, pdb_ in enumerate(sorted_pdbs):
            ag_ = prd.parsePDB(pdb_)
            ensemble.addCoordset(ag_.getCoords())
            os.remove(pdb_)
        ensemble.setAtoms(ag_)
        prd.writePDB(join(self.out_poses, 'poses_coords.pdb'), ensemble)

    def yield_filter_sort(self):

        # Get all scores
        scores_path = next(
            cmn.recursive_finder('poses_scores.txt', self.out_poses))
        with open(scores_path, 'rt') as inp:
            file_string = inp.readlines()
        scores = [x.strip() for x in file_string]

        scores_path_out = join(self.out_poses, 'scores.txt')
        with open(scores_path_out, 'wt') as out:
            for i, x in enumerate(scores):
                out.write(f'{i}    {x}\n')

        # Get all poses
        pdb_path = next(
            cmn.recursive_finder('poses_coords.pdb', self.out_poses))
        ensemble = root.Molecule(pdb_path).get_ensemble()
        out_name = join(self.out_poses, 'poses.pdb')
        prd.writePDB(out_name, ensemble)

        # Get filtered indices
        rec_kdt = ckd(self.rec.getCoords())
        lig_parsed = root.Molecule(pdb_path).parse()
        filtered_indices, filtered_ligs = root.get_filtered_indices(rec_kdt,
                                                                    lig_parsed)

        # Get filtered poses
        lig_filtered = prd.Ensemble()
        [lig_filtered.addCoordset(x.getCoords()) for x in filtered_ligs]
        lig_filtered.setAtoms(filtered_ligs[0])
        prd.writePDB(join(self.out_poses, 'poses_filtered.pdb'), lig_filtered)

        # Get filtered scores
        with open(join(self.out_poses, 'scores_filtered.txt'), 'wt') as out:
            for index in filtered_indices:
                out.write(f'{index}    {scores[index]}\n')

    def runner(self):
        self.build_hierarchy()
        self.prepare_ligand()
        self.prepare_receptor()
        self.prepare_gpf()
        self.compute_maps()
        self.prepare_dpf()
        self.dock()
        self.output_poses()
        self.convert_poses()
        self.create_pdbqt()
        self.output_scores()
        self.yield_filter_sort()

    def get_commands(self):
        return None, None


def parse_arguments():
    # Initializing argparse ---------------------------------------------------
    desc = '\nWrapper for ad4'
    parser = argparse.ArgumentParser(prog='ad4',
                                     description=desc,
                                     add_help=True,
                                     epilog='As simple as that ;)',
                                     allow_abbrev=False)
    # Arguments: loading trajectory -------------------------------------------
    all_args = parser.add_argument_group(title='Trajectory options')

    all_args.add_argument('-ad4_path', dest='ad4', action='store',
                          help='Path to ad4 root dir', type=str,
                          metavar='ad4', required=True)

    all_args.add_argument('-adt_path', dest='adt', action='store',
                          help='Path to adt root dir', type=str,
                          required=False, metavar='adt', default=None)

    all_args.add_argument('-babel_path', dest='babel', action='store',
                          help='Path to babel root dir', type=str,
                          required=True,
                          default=0, metavar='babel')

    all_args.add_argument('-rec_pdb', dest='rec_pdb', action='store',
                          help='Path to the rec.pdb', type=str, required=True,
                          default=None, metavar='rec_pdb')

    all_args.add_argument('-lig_pdb', dest='lig_pdb', action='store',
                          help='Path to the lig.pdb', type=str, required=True,
                          default=None, metavar='rec_pdb')

    all_args.add_argument('-n_poses', dest='n_poses', action='store',
                          help='Number of requested poses', type=int,
                          required=True, default=None, metavar='stride')

    all_args.add_argument('-rmsd_tol', dest='rmsd_tol', action='store',
                          help='rmsd cutoff for clustering', type=float,
                          required=True, default=None, metavar='rmsd_tol')

    all_args.add_argument('-exh', dest='exh', action='store',
                          help='exhaustiveness', type=float,
                          required=True, default=None, metavar='exh')

    all_args.add_argument('-odir', action='store', dest='odir',
                          help='Output directory to store analysis',
                          type=str, required=True, default=None,
                          metavar='odir')
    user_inputs = parser.parse_args()
    return user_inputs


def main():
    args = parse_arguments()  # for calling as cli

    ad4_obj = AutoDock4(args.ad4, args.adt, args.babel,
                        args.ad4, args.rec_pdb, args.lig_pdb,
                        args.n_poses, args.rmsd_tol, [], args.exh,
                        args.odir)
    ad4_obj.runner()


if __name__ == '__main__':
    main()

# =============================================================================
#
# =============================================================================
