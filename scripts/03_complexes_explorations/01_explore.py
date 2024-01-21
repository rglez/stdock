# Created by roy.gonzalez-aleman at 26/12/2023
import os
import subprocess as sp
from os.path import join

from tqdm import tqdm

import commons as cmn
import root
from programs.gnina import Gnina
from programs.plants import Plants
from programs.qvinaw import QvinaW
from programs.smina import Smina
from programs.vina import Vina


# todo: ensure correct ordering will be used in rmsd calculations (AD4 issue)
# todo: correct autodock4 exhaustiveness
# todo: add ability to restart in run_commands


def run_plants(plants_exe, rec_mol2, lig_mol2, n_poses, rmsd_tol, out_dir):
    """
    Run docking calculations on selected case using PLANTS
    Args:
        plants_exe: path to the PLANTS executable
        rec_mol2: receptor in mol2 format provided by SPORES
        lig_mol2: ligand in mol2 format provided by SPORES
        n_poses: number of requested poses
        rmsd_tol: rmsd cutoff for poses internal clustering
        out_dir: output directory where to put PLANTS files
    """
    plants_odir = join(out_dir, 'plants')
    plants_scores = ['plp', 'plp95', 'chemplp']
    plants_levels = ['speed1', 'speed2', 'speed4']

    PlantsObj = Plants(plants_exe, rec_mol2, lig_mol2, n_poses, rmsd_tol,
                       plants_scores, plants_levels, plants_odir)
    PlantsObj.run_commands()
    PlantsObj.yield_filter_sort()


def run_vina(vina_exe, rec_qt, lig_qt, n_poses, rmsd_tol, out_dir):
    """
    Run docking calculations on selected case using vina
    Args:
        vina_exe: path to the VINA executable
        rec_qt: receptor in pdbqt format provided by AutoDock Tools
        lig_qt: ligand in pdbqt format provided by AutoDock Tools
        n_poses: number of requested poses
        rmsd_tol: rmsd cutoff for poses internal clustering
        out_dir: output directory where to put PLANTS files
    """
    vina_odir = join(out_dir, 'vina')
    vina_scores = ['vina', 'vinardo']
    vina_levels = [8, 80, 800]
    VinaObj = Vina(vina_exe, rec_qt, lig_qt, n_poses, rmsd_tol,
                   vina_scores, vina_levels, vina_odir)
    VinaObj.run_commands()
    VinaObj.yield_filter_sort()


def run_qvinaw(qvinaw_exe, rec_qt, lig_qt, n_poses, rmsd_tol, out_dir):
    qvinaw_odir = join(out_dir, 'qvinaw')
    qvinaw_scores = []
    qvinaw_levels = [8, 80, 800]
    QvinawObj = QvinaW(qvinaw_exe, rec_qt, lig_qt, n_poses, rmsd_tol,
                       qvinaw_scores, qvinaw_levels, qvinaw_odir)
    QvinawObj.run_commands()
    QvinawObj.yield_filter_sort()


def run_smina(smina_exe, rec_qt, lig_qt, n_poses, rmsd_tol, out_dir):
    smina_scores = ['vina', 'vinardo', 'dkoes_fast', 'dkoes_scoring']
    smina_levels = [8, 80, 800]
    SminaObj = Smina(smina_exe, rec_qt, lig_qt, n_poses, rmsd_tol,
                     smina_scores, smina_levels, out_dir)
    SminaObj.run_commands()
    SminaObj.yield_filter_sort()


def run_gnina_no_cnn(gnina_exe, rec_qt, lig_qt, n_poses, rmsd_tol, out_dir):
    gnina_odir = join(out_dir, 'gnina')
    gnina_scores = ['ad4_scoring', 'default', 'dkoes_fast', 'dkoes_scoring',
                    'dkoes_scoring_old', 'vina', 'vinardo']
    gnina_levels = [8, 80, 800]
    GninaObj = Gnina(gnina_exe, rec_qt, lig_qt, n_poses, rmsd_tol,
                     gnina_scores, gnina_levels, gnina_odir)
    GninaObj.run_commands()
    GninaObj.yield_filter_sort()


def run_ad4(python_exe, ad4_exe, ad4_path, adt_path, babel_path, lig_pdb,
            n_poses, rmsd_tol, out_dir):

    n_poses = 2000 if n_poses > 2000 else n_poses
    ad4_levels = {25000: 'low', 250000: 'medium', 2500000: 'high'}
    bench = ['/usr/bin/time', '-v']
    for exh in tqdm(ad4_levels, total=len(ad4_levels)):
        cmd = (
            f'{python_exe} {ad4_exe} -ad4_path {ad4_path} -adt_path {adt_path}'
            f' -babel_path {babel_path} -rec_pdb {rec_pdb} -lig_pdb {lig_pdb}'
            f' -n_poses {n_poses} -rmsd_tol {rmsd_tol} -exh {exh} -odir {out_dir}')

        cmd_list = cmd.split()
        cmd_run = sp.Popen(bench + cmd_list, text=True, stdout=sp.PIPE,
                           stderr=sp.PIPE)
        output, errors = cmd_run.communicate()

        odir = join(out_dir, 'autodock4', f'autodock4_{ad4_levels[exh]}')
        log_name = join(odir, odir.split(os.sep)[-1] + '.log')
        root.write_string(output + errors, log_name)


# =============================================================================
# User-defined parameters
# =============================================================================
num_poses = 5000
rmsd_tol = 1.0
# =============================================================================

# ==== Prepare folders hierarchy
proj_dir = cmn.proj_dir
inputs_dir = cmn.get_abs(
    'scripts/02_complexes_preparation/01_prepared_with_adt')
cases = os.listdir(inputs_dir)
top_out_dir = cmn.get_abs('scripts/03_complexes_explorations/01_explored')
cmn.makedir_after_overwriting(top_out_dir)

# ==== Running programs for each case
for case in cases:
    # Get the files needed to run
    lig_qt = join(inputs_dir, case, f'{case}_ligand.mol2.pdbqt')
    lig_mol2 = join(inputs_dir, case, f'{case}_ligand.mol2')
    lig_pdb = join(inputs_dir, case, f'{case}_ligand.mol2.pdb')
    rec_qt = join(inputs_dir, case, f'{case}_protein.pdbqt')
    rec_pdb = join(inputs_dir, case, f'{case}_protein.pdb')
    rec_mol2 = join(inputs_dir, case, f'{case}_protein.mol2')

    # Build per-case hierarchy
    case_out_dir = join(top_out_dir, case)
    cmn.makedir_after_overwriting(case_out_dir)

    # =========================================================================
    # Run selected programs
    # =========================================================================

    # todo: WTF with mol2 charges with spores?
    # todo: be careful with mol2 having Mg
    # todo: make a converter from any to any using pybel. Convert to pdb with
    #       prody then back to mol2
    # print(f'Running PLANTS on {case}')
    # run_plants(cmn.plants_exe, rec_mol2, lig_mol2, num_poses, rmsd_tol,
    #            case_out_dir)

    # print(f'Running VINA on {case}')
    # run_vina(cmn.vina_exe, rec_qt, lig_qt, num_poses, rmsd_tol, case_out_dir)

    # print(f'Running QVINA-W on {case}')
    # run_qvinaw(cmn.qvinaw_exe, rec_qt, lig_qt, num_poses, rmsd_tol,
    #            case_out_dir)

    # print(f'Running SMINA on {case}')
    # run_smina(cmn.smina_exe, rec_qt, lig_qt, num_poses, rmsd_tol,
    #           case_out_dir)

    # todo: gnina without cnn rescoring? both?
    # print(f'Running GNINA on {case}')
    # run_gnina_no_cnn(cmn.gnina_exe, rec_qt, lig_qt, num_poses, rmsd_tol,
    #                  case_out_dir)

    # todo: what about autodock4-gpu ?
    print(f'Running AUTODOCK4 on {case}')
    run_ad4(cmn.python_exe, cmn.ad4_exe, cmn.ad4_path, cmn.adtools_dir,
            cmn.babel_path, lig_pdb, num_poses, rmsd_tol, case_out_dir)

# %%
# =============================================================================
# DOCK6
# =============================================================================

# =============================================================================
# EDOCK
# =============================================================================
