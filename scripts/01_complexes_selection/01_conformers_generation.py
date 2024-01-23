# Created by roy.gonzalez-aleman at 22/12/2023
"""
Generates "n_conformers" for every ligand in "core_set" with a rmsd distance
of "rmsd_tol". In case of rigid ligands, it is expected fewer conformers to be
 generated"
"""

import os
from os.path import join, split
import commons as cmn
import proj_paths as pp
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom, rdMolAlign
from tqdm import tqdm


def conformer_generator(input_mol2, out_traj_pdb, rmsd_cut, num_conf):
    """
    Generates ligand conformers

    Args:
        input_mol2: input ligand in mol2 format
        out_traj_pdb: name of the output pdb-formatted conformers trajectory
        rmsd_cut: do not generate conformers under this rmsd distance
        num_conf: number of conformers to generate. In case of rigid ligands,
        it is expected fewer conformers than this number to be generated
    """
    # Add H atoms
    try:
        ref_molecule = Chem.AddHs(
            Chem.MolFromMol2File(input_mol2), addCoords=True)
    except:
        print(f'Problems parsing {input_mol2}')
        return None

    # Generate conformers
    param = Chem.rdDistGeom.ETKDGv3()
    param.pruneRmsThresh = rmsd_cut
    param.randomSeed = 910623
    conformers_ids = Chem.rdDistGeom.EmbedMultipleConfs(
        ref_molecule, num_conf, param)

    # Minimize conformers
    mp = AllChem.MMFFGetMoleculeProperties(
        ref_molecule, mmffVariant='MMFF94s')
    AllChem.MMFFOptimizeMoleculeConfs(
        ref_molecule, numThreads=0, mmffVariant='MMFF94s')

    # Output trajectory
    w = Chem.PDBWriter(out_traj_pdb)
    res = []
    for cid in conformers_ids:
        ff = Chem.AllChem.MMFFGetMoleculeForceField(
            ref_molecule, mp, confId=cid)
        e = ff.CalcEnergy()
        res.append((cid, e))

    # Sort by energy
    sorted_res = sorted(res, key=lambda x: x[1])
    Chem.rdMolAlign.AlignMolConformers(ref_molecule)
    for cid, e in sorted_res:
        ref_molecule.SetProp('CID', str(cid))
        ref_molecule.SetProp('Energy', str(e))
        w.write(ref_molecule, confId=cid)
    w.close()

    # Output topology (the input mol2 file)
    root, base = split(out_traj_pdb)
    ref_molecule = Chem.AddHs(
        Chem.MolFromMol2File(input_mol2), addCoords=True)
    w = Chem.PDBWriter(join(root, f'top_{base}'))
    w.write(ref_molecule)
    w.close()


# =============================================================================
# User-specified arguments
# =============================================================================
rmsd_tol = 0.5
n_conformers = 100
# =============================================================================

# ==== Prepare folders hierarchy
root_dir = pp.proj_dir
core_set = join(root_dir, 'data/external/coreset')
output_dir = join(
    root_dir, 'scripts/01_complexes_selection/01_conformers_generation')
cmn.makedir_after_overwriting(output_dir)

# ==== Generate conformers
target = '*_ligand.mol2'
ligands = list(cmn.recursive_finder(target, core_set))

for ligand in tqdm(ligands, total=len(ligands)):
    # Get label
    label = ligand.split(os.sep)[-2]
    out_traj_name = join(output_dir, f'{label}_conformers.pdb')
    # Generate conformers
    conformer_generator(ligand, out_traj_name, rmsd_tol, n_conformers)
