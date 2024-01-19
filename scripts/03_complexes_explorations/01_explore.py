# Created by roy.gonzalez-aleman at 26/12/2023
from os.path import join

from programs.plants import Plants


# todo: ensure correct ordering will be used in rmsd calculations (AD4 issue)
# todo: correct autodock4 exhaustiveness
# todo: add ability to restart in run_commands


def run_plants(plants_exe, spores_exe, rec_mol2, lig_mol2, n_poses, rmsd_tol, out_dir):
    plants_odir = join(out_dir, 'plants')
    plants_scores = ['plp', 'plp95', 'chemplp']
    plants_levels = ['speed1', 'speed2', 'speed4']

    PlantsObj = Plants(plants_exe, spores_exe, rec_mol2, lig_mol2, n_poses, rmsd_tol,
                       plants_scores, plants_levels, plants_odir)
    PlantsObj.run_commands()
    PlantsObj.yield_filter_sort()


# =============================================================================
# User-defined parameters
# =============================================================================
# num_poses = 5000
# rmsd_tol = 1.0
# =============================================================================

# ==== Prepare folders hierarchy
# proj_dir = cmn.proj_dir
# inputs_dir = join(
#     proj_dir, 'scripts/02_complexes_preparation/01_prepared_with_adt')
# cases = os.listdir(inputs_dir)
# top_out_dir = join(proj_dir, 'scripts/03_complexes_explorations/01_explored')
# shutil.rmtree(top_out_dir, ignore_errors=True)
# os.makedirs(top_out_dir)

# ==== Running programs for each case
# for case in cases:
#     Get the files needed to run
# lig_qt = join(inputs_dir, case, f'{case}_ligand.pdbqt')
# lig_mol2 = join(inputs_dir, case, f'{case}_ligand.mol2')
# rec_qt = join(inputs_dir, case, f'{case}_protein.pdbqt')
# rec_pdb = join(inputs_dir, case, f'{case}_protein.pdb')
#
# Build per-case hierarchy
# case_out_dir = join(top_out_dir, case)
# os.makedirs(case_out_dir)
#
# Run selected programs
# print(f'Running X on {case}')


# %%

# =============================================================================
# Vina
# =============================================================================
# from programs.vina import Vina
#
# vina_exe = "/home/roy.gonzalez-aleman/SoftWare/vina_1.2.5_linux_x86_64"
# vina_odir = join(case_out_dir, 'vina')
# vina_scores = ['vina', 'vinardo']
# vina_levels = [8, 80, 800]
# VinaObj = Vina(vina_exe, rec_qt, lig_qt, n_poses, rmsd_tol,
#                vina_scores, vina_levels, vina_odir)
# VinaCommands = VinaObj.get_commands()
# VinaObj.run_commands()
# VinaObj.yield_filter_sort()
# =============================================================================
# Qvina-W
# =============================================================================
# from programs.qvinaw import QvinaW
#
# qvinaw_exe = '/home/roy.gonzalez-aleman/SoftWare/qvina-w'
# qvinaw_odir = join(case_out_dir, 'qvinaw')
# qvinaw_scores = []
# qvinaw_levels = [8, 80, 800]
# QvinawObj = QvinaW(qvinaw_exe, rec_qt, lig_qt, n_poses, rmsd_tol,
#                    qvinaw_scores, qvinaw_levels, qvinaw_odir)
# QvinawCommands = QvinawObj.get_commands()
# QvinawObj.run_commands()
# QvinawObj.yield_filter_sort()
# =============================================================================
# Smina
# =============================================================================
# from programs.smina import Smina
#
# smina_exe = "/home/roy.gonzalez-aleman/SoftWare/smina.static"
# smina_odir = join(case_out_dir, 'smina')
# smina_scores = ['vina', 'vinardo', 'dkoes_fast', 'dkoes_scoring']
# smina_levels = [8, 80, 800]
# SminaObj = Smina(smina_exe, rec_qt, lig_qt, n_poses, rmsd_tol,
#              smina_scores, smina_levels, smina_odir)
# SminaObj.run_commands()
# SminaObj.yield_filter_sort()

# =============================================================================
# Autodock4
# =============================================================================
# import root
# from tqdm import tqdm
# import subprocess as sp
#
# ad4_path = '/home/roy.gonzalez-aleman/SoftWare/autodock/x86_64Linux2/'
# adt_path = '/home/roy.gonzalez-aleman/SoftWare/autodock/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/'
# babel_path = '/home/roy.gonzalez-aleman/miniconda3/bin/obabel'
# ad4_exe = '/home/roy.gonzalez-aleman/RoyHub/stdock/src/stdock/ad4.py'
# ad4_levels = {250: 'low', 2500: 'medium', 25000: 'high'}
#
# bench = ['/usr/bin/time', '-v']
# for exh in tqdm(ad4_levels, total=len(ad4_levels)):
#     cmd = (f'{py_exe} {ad4_exe} -ad4_path {ad4_path} -adt_path {adt_path}'
#            f' -babel_path {babel_path} -rec_pdb {rec_pdb} -lig_pdb {lig_pdb}'
#            f' -n_poses {n_poses} -rmsd_tol {rmsd_tol} -exh {exh} -odir {case_out_dir}')
#
#     cmd_list = cmd.split()
#     cmd_run = sp.Popen(bench + cmd_list, text=True, stdout=sp.PIPE,
#                        stderr=sp.PIPE)
#     output, errors = cmd_run.communicate()
#     odir = join(case_out_dir, 'autodock4', f'autodock4_{ad4_levels[exh]}')
#     log_name = join(odir, odir.split(os.sep)[-1] + '.log')
#     root.write_string(output + errors, log_name)
# =============================================================================
# DOCK6
# =============================================================================
