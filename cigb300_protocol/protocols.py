# Created by roy.gonzalez-aleman at 26/02/2024
import os
from os.path import split, join, basename

import prody as prd
from tqdm import tqdm

import commons as cmn

out_dir = '/home/roy.gonzalez-aleman/RoyHub/stdock/cigb300_protocol/04_minimize_complexes/04_minimize_complex/minimized/'
complex_dir = '/home/roy.gonzalez-aleman/RoyHub/stdock/cigb300_protocol/04_minimize_complexes/03_parametrize_complex/parametrized/'
rec_psf = '/home/roy.gonzalez-aleman/RoyHub/stdock/cigb300_protocol/01_minimize_receptor/ck2_alpha.psf'
rec_pdb = '/home/roy.gonzalez-aleman/RoyHub/stdock/cigb300_protocol/03_docking_lig_conformers/stdock-CIGB300-389/DOCKING/lightdock_ck2_alpha_minim.pdb'

rec_parsed = prd.parsePDB(rec_pdb)
rec_parsed.setSegnames('R')
prd.writePDB(rec_pdb, rec_parsed)


def write_minimize_cmpx_tcl(complex_psf, complex_pdb, output_dir):
    dir_name, base_name_raw = split(complex_pdb)
    base_name = base_name_raw.split('.')[0]

    # ----
    script = f"""
structure         {complex_psf}
coordinates       {complex_pdb}

set outputname     complex
firsttimestep      0
paraTypeCharmm     on


parameters     /home/roy.gonzalez-aleman/RoyHub/stdock/cigb300_protocol/01_minimize_receptor/01_parametrize_receptor/toppar/toppar22/par_all22_prot.inp
parameters     /home/roy.gonzalez-aleman/RoyHub/stdock/cigb300_protocol/01_minimize_receptor/01_parametrize_receptor/toppar/toppar36_mod/par_all36_cgenff.prm
parameters     /home/roy.gonzalez-aleman/RoyHub/stdock/cigb300_protocol/01_minimize_receptor/01_parametrize_receptor/toppar/toppar36_mod/par_all36m_prot.prm
parameters     /home/roy.gonzalez-aleman/RoyHub/stdock/cigb300_protocol/01_minimize_receptor/01_parametrize_receptor/toppar/toppar36_mod/toppar_all36_prot_modify_res.str
parameters     /home/roy.gonzalez-aleman/RoyHub/stdock/cigb300_protocol/01_minimize_receptor/01_parametrize_receptor/toppar/toppar36_mod/par_all36_na.prm
parameters     /home/roy.gonzalez-aleman/RoyHub/stdock/cigb300_protocol/01_minimize_receptor/01_parametrize_receptor/toppar/toppar36_mod/par_all36_carb.prm
parameters     /home/roy.gonzalez-aleman/RoyHub/stdock/cigb300_protocol/01_minimize_receptor/01_parametrize_receptor/toppar/toppar36_mod/par_all36_lipid.prm


# force field ----------------------------------------------------
exclude            scaled1-4
1-4scaling         1.0
cutoff             12.
switching          on
switchdist         10.
pairlistdist       13.5
# integrator -----------------------------------------------------
timestep            1.0
rigidBonds          water
nonbondedFreq       1
fullElectFrequency  2
stepspercycle       20
# temperature ----------------------------------------------------
temperature         0
# output ---------------------------------------------------------
outputName          $outputname
restartfreq         1000
dcdfreq             50
xstFreq             1000
outputEnergies      1000
outputPressure      1000
outputTiming       1000
# run ------------------------------------------------------------
minimize            5000
"""
    # ----
    script_name = join(output_dir, f'minimize_complex_{base_name}.conf')
    with open(script_name, 'wt') as sn:
        sn.write(script)
    return script_name


ligs_psf_raw = list(cmn.recursive_finder('complex*psf', complex_dir))
ligs_pdb_raw = list(cmn.recursive_finder('complex*pdb', complex_dir))
assert len(ligs_psf_raw) == len(ligs_pdb_raw)

sort_files = lambda x: int(basename(x).split('.')[0].split('_')[-1])
ligs_psf = sorted(ligs_psf_raw, key=sort_files)
ligs_pdb = sorted(ligs_pdb_raw, key=sort_files)
zipped = zip(ligs_psf, ligs_pdb)

command = '/home/roy.gonzalez-aleman/SoftWare/NAMD_2.14_Linux-x86_64-multicore/namd2 +p4 {} > /dev/null'

# todo: uncomment for ab initio jobs
# minimize_commands = []
# for complex_psf, complex_pdb in zipped:
#     index = basename(complex_psf).split('.')[0].split('_')[-1]
#     odir = join(out_dir, index)
    # cmn.makedir_after_overwriting(odir)
    # script_path = write_minimize_cmpx_tcl(complex_psf, complex_pdb, odir)
    # cmd = command.format(script_path)
    # minimize_commands.append(cmd)



for task in tqdm(minimize_commands[487:], total=len(minimize_commands) - 487):
    os.system(task)
