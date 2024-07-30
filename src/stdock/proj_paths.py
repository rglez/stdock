# Created by rglez at 22/01/2024
from os.path import join

import stdock.commons as cmn

# ==== Path to the project root dir ===========================================
proj_dir = '/home/gonzalezroy/RoyHub/stdock'
cmn.check_path(proj_dir)

# ==== Project paths ==========================================================
plants_exe = join(proj_dir, 'programs/PLANTS1.2_64bit')
cmn.check_path(plants_exe)

spores_exe = join(proj_dir, 'programs/SPORES_64bit')
cmn.check_path(spores_exe)

vina_exe = join(proj_dir, 'programs/vina_1.2.5_linux_x86_64')
cmn.check_path(vina_exe)

qvinaw_exe = join(proj_dir, 'programs/qvina-w')
cmn.check_path(qvinaw_exe)

smina_exe = join(proj_dir, 'programs/smina.static')
cmn.check_path(smina_exe)

gnina_exe = join(proj_dir, 'programs/gnina')
cmn.check_path(gnina_exe)

ad4_exe = join(proj_dir, 'programs/ad4.py')
cmn.check_path(ad4_exe)

# ==== AutoDock4 paths ========================================================
ad4_root = '/home/gonzalezroy/SoftWare/autodock/'
cmn.check_path(ad4_root)

ad4_path = join(ad4_root, 'x86_64Linux2/')
cmn.check_path(ad4_path)

pythonsh_exe = join(ad4_root, 'mgltools_x86_64Linux2_1.5.7/bin/pythonsh')
cmn.check_path(pythonsh_exe)

adtools_dir = join(ad4_root,
                   'mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/')
cmn.check_path(adtools_dir)

# ==== Other paths ============================================================
python_exe = '/home/gonzalezroy/anaconda3/envs/stdock/bin/python'
cmn.check_path(python_exe)

babel_path = '/usr/bin/obabel'
cmn.check_path(babel_path)
