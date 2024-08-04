# Created by roy.gonzalez-aleman at 18/11/2023
"""
Compute std_score
"""
from os.path import basename

import numpy as np
import numpy_indexed as npi
import prody as prd
from scipy.spatial import cKDTree as ckd

import stdock.commons as cmn

epitope_str = 'EPITOPE-MAPPING*.txt'


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


def get_aliphatic_hydrogens(parsed_ag):
    # Get all bonds
    if not parsed_ag.getBonds():
        bonds = parsed_ag.inferBonds(set_bonds=False)
    else:
        bonds = parsed_ag.getBonds()

    # Get indices of H bonded to C
    aliphatic = []
    for at1, at2 in bonds:
        merged = f'{at1.getName()[0]}{at2.getName()[0]}'
        if merged == 'CH':
            aliphatic.append(at2.getIndex())
        elif merged == 'HC':
            aliphatic.append(at1.getIndex())
    return aliphatic


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


class STDScorer:
    """
    Score poses using std_score
    """

    def __init__(self, std_epitopes, topology, trajectory, multiplicity='1RXL',
                 max_dist=6, chain_rec='R', **kwargs):

        # Parse arguments
        self.epitopes = std_epitopes
        self.max_dist = max_dist
        self.topology = cmn.check_path(topology)
        self.trajectory = cmn.check_path(trajectory)
        self.chain_rec = chain_rec

        # Parse multiplicity-dependent arguments
        self.multiplicity = check_multiplicity(multiplicity)
        if self.multiplicity == '1RXL':
            self.rec_path = cmn.check_path(kwargs['rec_path'])
            rec = cmn.Molecule(self.rec_path).parse()[0]
            hc = get_aliphatic_hydrogens(rec)
            rec_hag = rec.select(f"index {' '.join([str(x) for x in hc])}")
            self.rec_htree = ckd(rec_hag.getCoords())
        self.qt_idxs = kwargs.get('qt_indices', None)

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

            # Re-naming epitope mapping taking into account pdbqt reordering
            matrix_labels = self.epitopes.sorted_names[lig_conc]
            matrix_idxs = [int(x.split('_')[1]) - 1 for x in matrix_labels]
            matrix_serial_int = npi.indices(self.qt_idxs, matrix_idxs) + 1
            matrix_serial_str = [str(x) for x in matrix_serial_int]
            matrix_labels = [f'{matrix_labels[i].split("_")[0]}_{x}'
                             for i, x in enumerate(matrix_serial_str)]

            # Parse topology & trajectory
            topo = prd.parsePDB(self.topology)
            traj = prd.Trajectory(self.trajectory)
            sele_epitope = topo.select(f"serial {' '.join(matrix_serial_str)}")
            sele_rec = topo.select(f'chain {self.chain_rec}')
            traj.setAtoms(sele_epitope)
            if not all([x.startswith('H') for x in sele_epitope.getNames()]):
                raise ValueError("Problems selecting the ligand's hydrogen")

            lig_atoms = traj.getAtoms()
            lig_names_all = lig_atoms.getNames()
            lig_serials = [str(x) for x in lig_atoms.getSerials()]
            lig_names = ['_'.join(x) for x in zip(lig_names_all, lig_serials)]
            lig_names = np.asarray(lig_names)

            # Gather std_score matrix info
            matrix = self.epitopes.score_matrices[lig_conc]
            rows = np.arange(matrix.shape[0])
            cols_inverse = rows[::-1]
            max_score = matrix[cols_inverse, rows].sum() / 2

            # Compute scores
            scores = []
            traj.reset()
            for frame in traj:
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
# input_dir = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/00-CASE_STUDY/HuR/M9/stdock-M9/'
# epitopes = STDEpitope(input_dir)
# topology = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/00-CASE_STUDY/HuR/M9/stdock-M9/topology.pdb'
# trajectory = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/00-CASE_STUDY/HuR/M9/stdock-M9/trajectory.dcd'
# rec_path = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/00-CASE_STUDY/HuR/M9/stdock-M9/receptor.pdb'
# self = STDScorer(epitopes, topology, trajectory, rec_path=rec_path)

# %% XRXL case
# input_dir = '/home/roy.gonzalez-aleman/RoyHub/stdock/cigb300_protocol/04_minimize_complexes/05_stdscore/'
# epitopes = STDEpitope(input_dir)
# topology = '/home/roy.gonzalez-aleman/RoyHub/stdock/cigb300_protocol/04_minimize_complexes/03_parametrize_complex/parametrized/complex_frame_0.pdb'
# trajectory = '/home/roy.gonzalez-aleman/RoyHub/stdock/cigb300_protocol/04_minimize_complexes/04_minimize_complex/minimized/0/complex.dcd'
# self = STDScorer(epitopes, topology, trajectory, multiplicity='XRXL',
#                  chain_lig='A', max_dist=7)

# Workflow

# input_dir = '/home/roy.gonzalez-aleman/RoyHub/stdock/cigb300_protocol/04_minimize_complexes/05_stdscore/'
# epitopes = STDEpitope(input_dir)
# topo_dir = '/home/roy.gonzalez-aleman/RoyHub/stdock/cigb300_protocol/04_minimize_complexes/03_parametrize_complex/parametrized/'
# traj_dir = '/home/roy.gonzalez-aleman/RoyHub/stdock/cigb300_protocol/04_minimize_complexes/04_minimize_complex/minimized'
#
# topologies_raw = list(cmn.recursive_finder('complex_frame*pdb', topo_dir))
# sort_topologies = lambda x: int(basename(x).split('.')[0].split('_')[-1])
# topologies = sorted(topologies_raw, key=sort_topologies)
#
# trajectories_raw = list(cmn.recursive_finder('complex.dcd', traj_dir))
# sort_trajs = lambda x: int(x.split(os.sep)[-2])
# trajectories = sorted(trajectories_raw, key=sort_trajs)
#
# for i, traj in enumerate(trajectories):
#     self = STDScorer(epitopes,
#                      topologies[i],
#                      traj,
#                      multiplicity='XRXL',
#                      chain_lig='A',
#                      max_dist=7)
#     last_score = self.scores[0.5][-1]
#     if last_score < 0.5:
#         print(traj, self.scores[0.5])
#
