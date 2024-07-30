# Created by roy.gonzalez-aleman at 26/12/2023
import os
import subprocess as sp
from os.path import join

from tqdm import tqdm

import commons as cmn
import proj_paths as pp
from programs.gnina import Gnina
from programs.plants import Plants
from programs.qvinaw import QvinaW
from programs.smina import Smina
from programs.vina import Vina


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
    vina_levels = [8, 64, 256]
    VinaObj = Vina(vina_exe, rec_qt, lig_qt, n_poses, rmsd_tol,
                   vina_scores, vina_levels, vina_odir)
    VinaObj.run_commands()
    VinaObj.yield_filter_sort()


def run_qvinaw(qvinaw_exe, rec_qt, lig_qt, n_poses, rmsd_tol, out_dir):
    qvinaw_odir = join(out_dir, 'qvinaw')
    qvinaw_scores = []
    qvinaw_levels = [8, 64, 256]
    QvinawObj = QvinaW(qvinaw_exe, rec_qt, lig_qt, n_poses, rmsd_tol,
                       qvinaw_scores, qvinaw_levels, qvinaw_odir)
    QvinawObj.run_commands()
    QvinawObj.yield_filter_sort()


def run_smina(smina_exe, rec_qt, lig_qt, n_poses, rmsd_tol, out_dir):
    smina_scores = ['vina', 'vinardo', 'dkoes_fast', 'dkoes_scoring']
    smina_levels = [8, 64, 256]
    SminaObj = Smina(smina_exe, rec_qt, lig_qt, n_poses, rmsd_tol,
                     smina_scores, smina_levels, out_dir)
    SminaObj.run_commands()
    SminaObj.yield_filter_sort()


def run_gnina_no_cnn(gnina_exe, rec_qt, lig_qt, n_poses, rmsd_tol, out_dir):
    gnina_odir = join(out_dir, 'gnina')
    gnina_scores = ['ad4_scoring', 'default']
    # 'dkoes_fast', 'dkoes_scoring', 'dkoes_scoring_old', 'vina', 'vinardo'
    # should all give the same results as smina
    gnina_levels = [8, 64, 256]
    GninaObj = Gnina(gnina_exe, rec_qt, lig_qt, n_poses, rmsd_tol,
                     gnina_scores, gnina_levels, gnina_odir)
    GninaObj.run_commands()
    GninaObj.yield_filter_sort()


def run_ad4(python_exe, ad4_exe, ad4_path, adt_path, babel_path, lig_pdb,
            n_poses, rmsd_tol, out_dir):
    n_poses = 2000 if n_poses > 2000 else n_poses
    # ad4_levels = {125000: 'low', 250000: 'medium', 500000: 'high'}
    ad4_levels = {500000: 'high'}
    bench = ['/usr/bin/time', '-v']
    for exh in tqdm(ad4_levels, total=len(ad4_levels)):
        cmd = (
            f'{python_exe} {ad4_exe} -ad4_path {ad4_path} -adt_path {adt_path}'
            f' -babel_path {babel_path} -rec_pdb {rec_pdb} -lig_pdb {lig_pdb}'
            f' -n_poses {n_poses} -rmsd_tol {rmsd_tol} -exh {exh} -odir {out_dir}')

        cmd_list = ' '.join(bench) + f' {cmd}'
        cmd_run = sp.Popen(cmd_list, text=True, stdout=sp.PIPE,
                           stderr=sp.PIPE, shell=True)
        output, errors = cmd_run.communicate()

        odir = join(out_dir, 'autodock4', f'autodock4_{ad4_levels[exh]}')
        log_name = join(odir, odir.split(os.sep)[-1] + '.log')
        cmn.write_string(output + errors, log_name)


def update_pickle(run_state_pickle, case=None, program=None):
    if os.path.exists(run_state_pickle):
        data = cmn.unpickle_from_file(run_state_pickle)
        if case and program:
            data[case][program] = True
        cmn.pickle_to_file(data, run_state_pickle)
    else:
        data = cmn.recursive_defaultdict()
        cmn.pickle_to_file(data, run_state_pickle)
    return data


# %%===========================================================================
# User-defined parameters
# =============================================================================
num_poses = 10000
rmsd_tol = 1.0
# =============================================================================

# Prepare folders hierarchy
proj_dir = pp.proj_dir
inputs_dir = join(proj_dir, 'scripts/02_complexes_preparation/01_prepared')
outputs_dir = join(proj_dir, 'scripts/03_complexes_explorations/01_explored')
if not os.path.exists(outputs_dir):
    os.makedirs(outputs_dir)
else:
    print('Restarting computations')

# Create a running state log for restarting jobs
run_state_pickle = join(outputs_dir, 'run_state.pickle')
run_state = update_pickle(run_state_pickle)

# Running programs for each case
# for case in (cases := os.listdir(inputs_dir)):
for case in ['2vkm', '4tmn']:
    # Build per-case hierarchy
    case_out_dir = join(outputs_dir, case)
    if not os.path.exists(case_out_dir):
        os.makedirs(case_out_dir)
    else:
        print(f'Restarting computations inside {case}')

    # Gather input files
    lig_pdb = join(inputs_dir, case, f'{case}_ligand.pdb')
    rec_pdb = join(inputs_dir, case, f'{case}_protein.pdb')
    lig_qt = join(inputs_dir, case, f'{case}_ligand.pdbqt')
    rec_qt = join(inputs_dir, case, f'{case}_protein.pdbqt')
    lig_spored = join(inputs_dir, case, f'{case}_ligand_spored.mol2')
    rec_spored = join(inputs_dir, case, f'{case}_protein_spored.mol2')

    # # Run PLANTS
    # if not run_state[case]['PLANTS']:
    #     print(f'Running PLANTS on {case}')
    #     run_plants(
    #         pp.plants_exe, rec_spored, lig_spored, num_poses, rmsd_tol,
    #         case_out_dir)
    #     run_state = update_pickle(run_state_pickle, case=case,
    #                               program='PLANTS')
    # else:
    #     print(f'Skipping PLANTS on {case} as already computed')

    # # Run VINA
    # if not run_state[case]['VINA']:
    #     print(f'Running VINA on {case}')
    #     run_vina(pp.vina_exe, rec_qt, lig_qt, num_poses, rmsd_tol,
    #              case_out_dir)
    #     run_state = update_pickle(run_state_pickle, case=case, program='VINA')
    # else:
    #     print(f'Skipping VINA on {case} as already computed')
    #
    # # Run QVINA-W
    # if not run_state[case]['QVINAW']:
    #     print(f'Running QVINAW on {case}')
    #     run_qvinaw(pp.qvinaw_exe, rec_qt, lig_qt, num_poses, rmsd_tol,
    #                case_out_dir)
    #     run_state = update_pickle(run_state_pickle, case=case,
    #                               program='QVINAW')
    # else:
    #     print(f'Skipping QVINAW on {case} as already computed')
    #
    # # Run SMINA
    # if not run_state[case]['SMINA']:
    #     print(f'Running SMINA on {case}')
    #     run_smina(pp.smina_exe, rec_qt, lig_qt, num_poses, rmsd_tol,
    #               case_out_dir)
    #     run_state = update_pickle(run_state_pickle, case=case, program='SMINA')
    # else:
    #     print(f'Skipping SMINA on {case} as already computed')
    #
    # # ==== Run GNINA
    # # todo: gnina without cnn rescoring? both?
    # if not run_state[case]['GNINA']:
    #     print(f'Running GNINA on {case}')
    #     run_gnina_no_cnn(pp.gnina_exe, rec_qt, lig_qt, num_poses, rmsd_tol,
    #                      case_out_dir)
    #     run_state = update_pickle(run_state_pickle, case=case, program='GNINA')
    # else:
    #     print(f'Skipping GNINA on {case} as already computed')
    #
    # # ==== Run AUTODOCK4
    # todo: what about autodock4-gpu ?
    if not run_state[case]['AD4']:
        print(f'Running AD4 on {case}')
        run_ad4(pp.python_exe, pp.ad4_exe, pp.ad4_path, pp.adtools_dir,
                pp.babel_path, lig_pdb, num_poses, rmsd_tol, case_out_dir)
        run_state = update_pickle(run_state_pickle, case=case, program='AD4')
    else:
        print(f'Skipping AD4 on {case} as already computed')
# %%
# =============================================================================
# DOCK6
# =============================================================================

# =============================================================================
# EDOCK
# =============================================================================
