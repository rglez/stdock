# Created by roy.gonzalez-aleman at 26/12/2023
from os.path import abspath, join

import docking

# =============================================================================
# General Parameters
# =============================================================================
root = abspath('.')
lig_qt = join(root,
              'scripts/dude_explorations/input_files/p97ND1/ligand.pdbqt')
lig_mol2 = join(root,
                'scripts/dude_explorations/input_files/p97ND1/ligand.mol2')
rec_qt = join(root,
              'scripts/dude_explorations/input_files/p97ND1/receptor.pdbqt')
rec_mol2 = join(root,
                'scripts/dude_explorations/input_files/p97ND1/receptor.mol2')

lig_pdb = '/home/roy.gonzalez-aleman/Tutorials/Bioinformatica/Posgrado_Bioinf/lig_Strytilcysteine.pdb'
rec_pdb = '/home/roy.gonzalez-aleman/Tutorials/Bioinformatica/Posgrado_Bioinf/rec_Eg5.pdb'

out_dir = join(root, 'scripts/dude_explorations/explorations/')
py_exe = '/home/roy.gonzalez-aleman/miniconda3/envs/stdock/bin/python'
n_poses = 2000
rmsd_tol = 1.0

# %%
# =============================================================================
# EDOCK
# =============================================================================
# =============================================================================
# DOCK6
# =============================================================================
# =============================================================================
# DOCK38
# =============================================================================
# =============================================================================
# LeDOCK
# =============================================================================
# =============================================================================
# Qvina-W
# =============================================================================
qvinaw_exe = '/home/roy.gonzalez-aleman/SoftWare/qvina-w'
qvinaw_odir = join(out_dir, 'qvinaw')
qvinaw_scores = []
qvinaw_levels = [8, 80, 800]
QvinawObj = docking.QvinaW(qvinaw_exe, rec_qt, lig_qt, n_poses, rmsd_tol,
                           qvinaw_scores, qvinaw_levels, qvinaw_odir)
QvinawCommands = QvinawObj.get_commands()

# =============================================================================
# Gnina
# =============================================================================
gnina_exe = '/home/roy.gonzalez-aleman/SoftWare/gnina'
gnina_odir = join(out_dir, 'gnina')
gnina_scores = ['ad4_scoring', 'default', 'dkoes_fast', 'dkoes_scoring',
                'dkoes_scoring_old', 'vina', 'vinardo']
gnina_levels = [8, 80, 800]
GninaObj = docking.Gnina(gnina_exe, rec_qt, lig_qt, n_poses, rmsd_tol,
                         gnina_scores, gnina_levels, gnina_odir)
GninaCommands = GninaObj.get_commands()

# =============================================================================
# Autodock4
# =============================================================================
ad4_exe = '/home/roy.gonzalez-aleman/RoyHub/stdock/src/stdock/ad4.py'
ad4_odir = join(out_dir, 'ad4')
ad4_levels = [2500000, 25000000, 250000000]
AutoDock4Commands = []
for level in ad4_levels:
    cmd = (
        f'{py_exe} {ad4_exe} {rec_pdb} {lig_pdb} {ad4_odir} {level} {n_poses}'
        f' {rmsd_tol}')
    AutoDock4Commands.append(cmd)

# =============================================================================
# PLANTS
# =============================================================================
plants_exe = '/home/roy.gonzalez-aleman/SoftWare/PLANTS1.2_64bit'
plants_odir = join(out_dir, 'plants')
plants_scores = ['plp', 'plp95', 'chemplp']
plants_levels = ['speed1', 'speed2', 'speed4']
PlantsObj = docking.Plants(plants_exe, rec_mol2, lig_mol2, n_poses, rmsd_tol,
                           plants_scores, plants_levels, plants_odir)
PlantsCommands = PlantsObj.get_commands()

# =============================================================================
# Vina
# =============================================================================
vina_exe = "/home/roy.gonzalez-aleman/SoftWare/vina_1.2.5_linux_x86_64"
vina_odir = join(out_dir, 'vina')
vina_scores = ['vina', 'vinardo']
vina_levels = [8, 80, 800]
VinaObj = docking.Vina(vina_exe, rec_qt, lig_qt, n_poses, rmsd_tol,
                       vina_scores, vina_levels, vina_odir)
VinaCommands = VinaObj.get_commands()

# =============================================================================
# Smina
# =============================================================================
smina_exe = "/home/roy.gonzalez-aleman/SoftWare/smina.static"
smina_odir = join(out_dir, 'smina')
smina_scores = ['vina', 'vinardo', 'dkoes_fast', 'dkoes_scoring']
smina_levels = [8, 80, 800]
SminaObj = docking.Smina(smina_exe, rec_qt, lig_qt, n_poses, rmsd_tol,
                         smina_scores, smina_levels, smina_odir)
SminaCommands = SminaObj.get_commands()
