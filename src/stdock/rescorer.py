# Created by roy.gonzalez-aleman at 18/11/2023
"""
Compute std_score
"""
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


class STDScorer:
    """
    Score poses using std_score
    """

    def __init__(self, lig_path, rec_path, poses_dcd, std_epitopes,
                 max_dist=6):

        # Parse arguments
        self.lig_path = cmn.check_path(lig_path)
        self.rec_path = cmn.check_path(rec_path)
        self.dcd_poses = cmn.check_path(poses_dcd)
        self.epitopes = std_epitopes
        self.max_dist = max_dist

        # Get scores
        rec_hag = cmn.Molecule(self.rec_path).parse()[0].hydrogen
        self.rec_htree = ckd(rec_hag.getCoords())
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
            serials = [x.split('_')[-1] for x in mapping]
            traj = prd.Trajectory(self.dcd_poses)
            topo = prd.parsePDB(self.lig_path)
            traj.setAtoms(topo.select(f"serial {' '.join(serials)}"))

            # Compute scores
            indices = [x - 1 for x in map(int, serials)]
            zipped = zip(topo.getNames()[indices], topo.getSerials()[indices])
            lig_names = np.asarray([f'{x[0]}_{x[1]}'for x in zipped])
            scores = []
            traj.reset()
            for frame in traj:
                lig_coords = frame.getCoords()
                distances, counterpart = self.rec_htree.query(lig_coords)
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


# %%
# receptor_path = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/TROLL8/p97ND1.pdbqt'
# poses_filtered = next(cmn.recursive_finder(poses_str, input_dir))
# poses = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/hpepdock/65cc7d529a1c3/hpepdock_all.pdb'
# receptor = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/hpepdock/65cc7d529a1c3/rec_65cc7d529a1c3.pdb'
# rec_path = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/stdock-CIGB300-lightdock2/lightdock_ck2_alpha.pdb'

# lig_path = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/stdock-CIGB300/DOCKING/lightdock_CIGB300_MINIM.pdb'
# rec_path = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/stdock-CIGB300/DOCKING/lightdock_ck2_alpha.pdb'
# poses_dcd = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/stdock-CIGB300/DOCKING/swarms.dcd'
#
# epitope_str = 'EPITOPE-MAPPING*.txt'
# input_dir = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/stdock-CIGB300_saved/'
# epitopes = STDEpitope(input_dir)
#
# self = STDScorer(lig_path, rec_path, poses_dcd, epitopes, max_dist=6)
