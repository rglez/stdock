# Created by roy.gonzalez-aleman at 18/11/2023
"""
Compute std_score
"""
import os
from os.path import basename

import numpy as np
import numpy_indexed as npi
import prody as prd
from scipy.spatial import cKDTree as ckd

import commons as cmn

epitope_str = 'EPITOPE-MAPPING*.txt'


class STDEpitope:
    """
    Parse STDock epitope mapping files & Compute the stdscore matrices
    """

    def __init__(self, input_dir):
        self.inp_dir = cmn.check_path(input_dir)

        # Get parsed mappings
        self.mappings_raw = list(cmn.recursive_finder(epitope_str, input_dir))
        self.parsed_mappings = self.parse_mappings()

        # Get score matrices
        self.sorted_names = self.get_sorted_names()
        self.score_matrices = self.get_score_matrices()

    def parse_mappings(self):
        """
        Parse epitope mappings inside the input dir

        Returns:
            mappings_dict: a dict of {lig_conc: std_mapping}
        """
        mappings_dict = {}
        for mapping in self.mappings_raw:
            lig_conc = float(basename(mapping).split('_')[1].split('.txt')[0])
            mappings_dict.update({lig_conc: {}})
            with open(mapping, 'rt') as mapp:
                for line in mapp:
                    splitted = line.split(':')
                    h_name, std_norm = splitted[0].strip(), float(splitted[1])
                    mappings_dict[lig_conc].update({h_name: std_norm})
        return mappings_dict

    def get_sorted_names(self):
        """
        Get ligand H names sorted from higher std normalizaed intensity for
        score matrix construction

        Returns:
            sorted_names: sorted ligand H names
        """
        sorted_names = {}
        for lig_conc in self.parsed_mappings:
            mapping = self.parsed_mappings[lig_conc]
            reordered = sorted(mapping, key=lambda x: mapping[x], reverse=True)
            sorted_names.update({lig_conc: np.asarray(reordered)})
        return sorted_names

    def get_score_matrices(self):
        """
        Get stds_core matrices for ligand poses recoring using std info

        Returns:
            score_matrices: std_score matrices
        """
        sorted_names = self.sorted_names
        score_matrices = {}
        for lig_conc in sorted_names:
            names = sorted_names[lig_conc]
            score_matrix = []
            for name in names:
                scores = []
                for alternate in names:
                    score = abs(self.parsed_mappings[lig_conc][name] -
                                self.parsed_mappings[lig_conc][alternate])
                    scores.append(score)
                score_matrix.append(scores)
            score_matrices.update({lig_conc: np.asarray(score_matrix)})
        return score_matrices


def check_multiplicity(multiplicity):
    """
    Check the multiplicity value provided is supported
    Args:
        multiplicity: 1RXL for a single receptor and multiple ligands
                      XRXL for as many receptors as ligands
    Returns:
        multiplicity: the passed value if it is supported else raise
    """
    if multiplicity not in (multiplicities := ['1RXL', 'XRXL']):
        raise ValueError(f'Multiplicity must be one of: {multiplicities}')
    return multiplicity


class STDScorer:
    """
    Score poses using std_score
    """

    def __init__(self, std_epitopes, topology, trajectory, multiplicity='1RXL',
                 max_dist=6, chain_rec='R', chain_lig='L', **kwargs):

        # Parse arguments
        self.epitopes = std_epitopes
        self.max_dist = max_dist
        self.topology = cmn.check_path(topology)
        self.trajectory = cmn.check_path(trajectory)
        self.chain_rec = chain_rec
        self.chain_lig = chain_lig

        # Parse multiplicity-dependent arguments
        self.multiplicity = check_multiplicity(multiplicity)
        if self.multiplicity == '1RXL':
            self.rec_path = cmn.check_path(kwargs['rec_path'])
            rec_hag = cmn.Molecule(self.rec_path).parse()[0].hydrogen
            self.rec_htree = ckd(rec_hag.getCoords())

        # Get scores
        self.scores = self.rescore_poses()

    def rescore_poses(self):
        """
        Re-score ligands poses using std_score matrices

        Returns:
            rescores: a dict of {lig_conc: std_score for each pose}
        """

        parsed_mappings = self.epitopes.parsed_mappings
        rescores = {}
        for lig_conc in parsed_mappings:
            rescores.update({lig_conc: {}})

            # Gather std_score matrix info
            score_matrices = self.epitopes.score_matrices
            matrix = score_matrices[lig_conc]
            rows = np.arange(matrix.shape[0])
            cols_inverse = rows[::-1]
            max_score = matrix[cols_inverse, rows].sum() / 2
            matrix_labels = self.epitopes.sorted_names[lig_conc]

            # Select only ligand atoms appearing on mappings
            mapping = parsed_mappings[lig_conc]
            serials_mapping_str = [x.split('_')[-1] for x in mapping]
            serials_mapping_int = [int(x) for x in serials_mapping_str]
            indices_mapping = [x - 1 for x in serials_mapping_int]

            topo = prd.parsePDB(self.topology)
            # serials_topo = topo.select(f'chain {self.chain_lig}').getSerials()
            # indices_topo = npi.indices(serials_topo, serials_mapping_int)
            # serials_real = serials_topo[indices_topo]
            traj = prd.Trajectory(self.trajectory)
            sele_epitope = topo.select(
                f"serial {' '.join(serials_mapping_str)}")
            sele_rec = topo.select(f'chain {self.chain_rec}')

            # Get H names correctly
            traj.setAtoms(sele_epitope)
            zipped = zip(topo.getNames()[indices_mapping],
                         topo.getSerials()[indices_mapping])
            lig_names = np.asarray([f'{x[0]}_{x[1]}' for x in zipped])
            if not all([x.startswith('H') for x in lig_names]):
                raise ValueError("Problems selecting the ligand's hydrogen")

            # Compute scores
            scores = []
            traj.reset()
            for frame in traj:
                traj.setAtoms(sele_epitope)
                lig_coords = frame.getCoords()

                if self.multiplicity == '1RXL':
                    distances, counterpart = self.rec_htree.query(lig_coords)
                else:
                    traj.setAtoms(sele_rec)
                    rec_htree = ckd(frame.getAtoms().hydrogen.getCoords())
                    distances, counterpart = rec_htree.query(lig_coords)

                if all(distances <= self.max_dist):
                    permutation = lig_names[distances.argsort()]
                    cols = npi.indices(matrix_labels, permutation)
                    score = matrix[cols, rows].sum() / 2 / max_score
                    scores.append(score)
                else:
                    scores.append(2)
            rescores.update({lig_conc: np.asarray(scores)})

        return rescores


# =============================================================================
# Debugging area
# =============================================================================

# %% 1RXL case
# rec_path = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/stdock-CIGB300/DOCKING/lightdock_ck2_alpha.pdb'
# input_dir = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/stdock-CIGB300-389'
# epitopes = STDEpitope(input_dir)
# topology = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/stdock-CIGB300-389/DOCKING/lightdock_cigb300_389.pdb'
# trajectory = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/stdock-CIGB300-389/DOCKING/swarms.dcd'
# self = STDScorer(epitopes, topology, trajectory, rec_path=rec_path)

# %% XRXL case
# input_dir = '/home/roy.gonzalez-aleman/RoyHub/stdock/cigb300_protocol/04_minimize_complexes/05_stdscore/'
# epitopes = STDEpitope(input_dir)
# topology = '/home/roy.gonzalez-aleman/RoyHub/stdock/cigb300_protocol/04_minimize_complexes/03_parametrize_complex/parametrized/complex_frame_0.pdb'
# trajectory = '/home/roy.gonzalez-aleman/RoyHub/stdock/cigb300_protocol/04_minimize_complexes/04_minimize_complex/minimized/0/complex.dcd'
# self = STDScorer(epitopes, topology, trajectory, multiplicity='XRXL',
#                  chain_lig='A', max_dist=7)

# Workflow

input_dir = '/home/roy.gonzalez-aleman/RoyHub/stdock/cigb300_protocol/04_minimize_complexes/05_stdscore/'
epitopes = STDEpitope(input_dir)
topo_dir = '/home/roy.gonzalez-aleman/RoyHub/stdock/cigb300_protocol/04_minimize_complexes/03_parametrize_complex/parametrized/'
traj_dir = '/home/roy.gonzalez-aleman/RoyHub/stdock/cigb300_protocol/04_minimize_complexes/04_minimize_complex/minimized'

topologies_raw = list(cmn.recursive_finder('complex_frame*pdb', topo_dir))
sort_topologies = lambda x: int(basename(x).split('.')[0].split('_')[-1])
topologies = sorted(topologies_raw, key=sort_topologies)

trajectories_raw = list(cmn.recursive_finder('complex.dcd', traj_dir))
sort_trajs = lambda x: int(x.split(os.sep)[-2])
trajectories = sorted(trajectories_raw, key=sort_trajs)

for i, traj in enumerate(trajectories):
    self = STDScorer(epitopes,
                     topologies[i],
                     traj,
                     multiplicity='XRXL',
                     chain_lig='A',
                     max_dist=7)
    last_score = self.scores[0.5][-1]
    if last_score < 0.5:
        print(traj, self.scores[0.5])

# =============================================================================
#
# =============================================================================
import matplotlib.pyplot as plt

cmn.generic_matplotlib(width=(10, 5))
traj_stds = ['''2.         2.         2.         2.         2.         2.
 2.         2.         2.         2.         2.         2.
 2.         2.         2.         2.         2.         2.
 2.         2.         2.         2.         2.         2.
 2.         2.         2.         2.         2.         2.
 2.         2.         2.         2.         2.         2.
 2.         2.         2.         2.         2.         2.
 2.         2.         2.         2.         2.         2.
 2.         2.         2.         2.         2.         2.
 2.         2.         2.         2.         2.         2.
 2.         0.39662447 0.39662447 0.39662447 0.39662447 0.39662447
 0.39662447 0.39662447 0.39662447 0.39662447 0.40084388 0.40084388
 0.40084388 0.40084388 0.40084388 0.40084388 0.40084388 0.40084388
 0.40084388 0.40084388 0.40084388 0.40084388 0.40084388 0.40084388
 0.40084388 0.40084388 0.40084388 0.40084388 0.40084388 0.40084388
 0.40084388 0.40084388 0.40084388 0.46413502 0.46413502 0.46413502
 0.40084388 0.40084388 0.40084388 0.40084388''',
             """0.72151899 0.80168776 0.78059072 0.78059072 0.78059072 0.70042194
              0.78059072 0.80168776 0.82278481 0.80168776 0.80168776 0.71729958
              0.64135021 0.59915612 0.59915612 0.69620253 0.73839662 0.82278481
              0.74261603 0.74261603 0.74261603 0.74261603 0.74261603 0.79746835
              0.73417722 0.65822785 0.57805907 0.53586498 0.47257384 0.47257384
              0.47257384 0.53586498 0.47257384 0.39662447 0.39662447 0.40084388
              0.40084388 0.40084388 0.37974684 0.35443038 0.45147679 0.45147679
              0.47257384 0.47257384 0.47257384 0.47257384 0.49367089 0.48945148
              0.56540084 0.48945148 0.41350211 0.41350211 0.41350211 0.38818565
              0.42194093 0.42194093 0.46413502 0.46413502 0.46413502 0.46413502
              0.46413502 0.46413502 0.46413502 0.38396624 0.31223629 0.31223629
              0.31223629 0.3164557  0.4556962  0.53586498 0.53586498 0.53586498
              0.53586498 0.53586498 0.39662447 0.2742616  0.2742616  0.2742616
              0.27004219 0.27004219 0.27004219 0.2742616  0.2742616  0.2742616
              0.2742616  0.2742616  0.2742616  0.2742616  0.27004219 0.2742616
              0.27004219 0.2742616  0.2742616  0.2742616  0.2742616  0.2742616
              0.2742616  0.2742616  0.2742616  0.2742616
             """,
             """0.4092827  0.38818565 0.4092827  0.4092827  0.43037975 0.51054852
              0.51054852 0.51054852 0.51054852 0.48945148 0.41772152 0.41772152
              0.41772152 0.41772152 0.49367089 0.41772152 0.41772152 0.43881857
              0.43881857 0.51476793 0.51476793 0.51476793 0.51898734 0.51898734
              0.43459916 0.43459916 0.43459916 0.43881857 0.43881857 0.43881857
              0.41772152 0.41772152 0.41772152 0.41350211 0.41350211 0.41350211
              0.41772152 0.41772152 0.41772152 0.42194093 0.42194093 0.42194093
              0.34177215 0.42194093 0.42194093 0.42194093 0.42194093 0.42194093
              0.34177215 0.42194093 0.42194093 0.42194093 0.42194093 0.42194093
              0.34177215 0.34177215 0.32067511 0.40084388 0.3164557  0.3164557
              0.32067511 0.32067511 0.32067511 0.32067511 0.32067511 0.32067511
              0.32067511 0.32067511 0.32067511 0.32067511 0.32067511 0.32067511
              0.32067511 0.32067511 0.32067511 0.32067511 0.32067511 0.32067511
              0.32067511 0.32067511 0.32067511 0.32067511 0.32067511 0.32067511
              0.32067511 0.3164557  0.32067511 0.32067511 0.32067511 0.32067511
              0.32067511 0.3164557  0.3164557  0.3164557  0.3164557  0.3164557
              0.3164557  0.32067511 0.3164557  0.32067511"""
             ]

y_points = []
for x in traj_stds:
    points = [float(x) for x in x.split()]
    y_points.append(points)
    plt.plot(points, lw=1, marker='.', ms=10)
plt.ylabel('STDScore', weight='bold')
plt.xlabel('Minim Step', weight='bold')
plt.savefig('std_trajs')
plt.close()
