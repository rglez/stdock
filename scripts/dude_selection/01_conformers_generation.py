# Created by roy.gonzalez-aleman at 22/12/2023
"""
This script is conceived to generate X conformers for each ligand in the DUDE
dataset.
"""
import os
from os.path import join, split

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdMolAlign

import commons as cmn


def conformer_generator(input_mol2, out_traj_pdb, rmsd_cut, num_conf):
    # Conformer generation
    try:
        ref_molecule = Chem.AddHs(Chem.MolFromMol2File(input_mol2), addCoords=True)
    except:
        print(f'Problems parsing {input_mol2}')
        return None

    param = rdDistGeom.ETKDGv3()
    param.pruneRmsThresh = rmsd_cut
    param.randomSeed = 910623
    conformers_ids = rdDistGeom.EmbedMultipleConfs(ref_molecule, num_conf,
                                                   param)

    # Conformer minimization
    mp = AllChem.MMFFGetMoleculeProperties(ref_molecule, mmffVariant='MMFF94s')
    AllChem.MMFFOptimizeMoleculeConfs(ref_molecule, numThreads=0,
                                      mmffVariant='MMFF94s')

    # Output trajectory
    w = Chem.PDBWriter(out_traj_pdb)
    res = []
    for cid in conformers_ids:
        ff = AllChem.MMFFGetMoleculeForceField(ref_molecule, mp, confId=cid)
        e = ff.CalcEnergy()
        res.append((cid, e))

    # Sorted by energy
    sorted_res = sorted(res, key=lambda x: x[1])
    rdMolAlign.AlignMolConformers(ref_molecule)
    for cid, e in sorted_res:
        ref_molecule.SetProp('CID', str(cid))
        ref_molecule.SetProp('Energy', str(e))
        w.write(ref_molecule, confId=cid)
    w.close()

    # Output topology
    root, base = split(out_traj_pdb)
    ref_molecule = Chem.AddHs(Chem.MolFromMol2File(input_mol2), addCoords=True)
    w = Chem.PDBWriter(join(root, f'top_{base}'))
    w.write(ref_molecule)
    w.close()


# =============================================================================
# User-specified arguments
# =============================================================================
dude_dir = '/home/roy.gonzalez-aleman/DataHub/data_raw/dude_all/'
target = 'crystal_ligand.mol2'

root_dir = os.path.abspath('.')
output_dir = join(root_dir, 'scripts/dude_selection/01_conformers_generation')
os.makedirs(output_dir, exist_ok=True)
rmsd_tol = 0.5
n_conformers = 100

# =============================================================================
# Conformers generation
# =============================================================================
ligands = list(cmn.recursive_finder(target, dude_dir))
for ligand in ligands:
    # Get label
    label = ligand.split(os.sep)[-2]
    out_traj_name = join(output_dir, f'{label}_conformers.pdb')

    # Generate conformers
    print(f'Generating {n_conformers} conformers with rmsd_tol of {rmsd_tol}'
          f' for target {label}')
    conformer_generator(ligand, out_traj_name, rmsd_tol, n_conformers)
