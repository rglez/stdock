# Created by roy.gonzalez-aleman at 18/11/2023
import docking
import numpy as np
import numpy_indexed as npi
import prody as prd
from matplotlib import pyplot as plt
from scipy.spatial import cKDTree as ckd


# todo: choose only aliphatic hydrogens


class STDScorer:
    def __init__(self, receptor, ligands, epitope_mapping):

        # Get the receptor kdtree
        self.rec_ag = receptor
        self.rec_hag = self.rec_ag.hydrogen
        self.rec_htree = ckd(self.rec_hag.getCoords())

        # Parse ligands and set a reference
        self.lig_ags = ligands
        self.ref_ligand = self.lig_ags[0]

        # Check all H atoms passed as epitope are in the lig molecular file
        self.raw_mapping = self.check_mapping(epitope_mapping)
        self.mapping = np.asarray([x for x in self.raw_mapping.keys()])

        # Select for scoring only H atoms declared in the mapping
        self.lig_hags = [ag.select(f'name {" ".join(self.mapping)}') for ag
                         in self.lig_ags]

        # Get the scoring matrix & the max score
        self.matrix_labels = self.get_reference()
        self.matrix = self.get_score_matrix()
        self.rows = np.arange(self.matrix.shape[0])
        self.max_score = self.get_max_score()

        # Score poses
        self.scores = self.score()

    def check_mapping(self, epitope_mapping):
        all_names = self.ref_ligand.getNames()
        for name in epitope_mapping:
            if name not in all_names:
                raise ValueError(
                    f'{name} declared in your mapping is not the passed pdb')
        return epitope_mapping

    def get_reference(self):
        ordered_names = sorted(
            self.raw_mapping, key=lambda x: int(self.raw_mapping[x]),
            reverse=True)
        return np.asarray(ordered_names)

    def get_score_matrix(self):
        score_matrix = []
        for name in self.matrix_labels:
            scores = []
            for alternate in self.matrix_labels:
                score = abs(mapping[name] - mapping[alternate])
                scores.append(score)
            score_matrix.append(scores)
        return np.asarray(score_matrix)

    def get_max_score(self):
        cols = npi.indices(self.matrix_labels, self.matrix_labels[::-1])
        return self.matrix[cols, self.rows].sum() / 2

    # def score_permutation(self, permutation):
    #     cols = npi.indices(self.mapping, permutation)
    #     return (self.matrix[cols, self.rows].sum() / 2) / self.max_score

    def score(self):
        lig_names = self.lig_hags[0].getNames()
        scores = []
        for lig in self.lig_hags:
            distances, counterpart = self.rec_htree.query(lig.getCoords())
            permutation = lig_names[distances.argsort()]

            cols = npi.indices(self.matrix_labels, permutation)

            score = self.matrix[cols, self.rows].sum() / 2 / self.max_score
            scores.append(score)
        return np.asarray(scores)


#
# =============================================================================
# Debugging area
# =============================================================================
# inputs
pdb_path = '/home/roy.gonzalez-aleman/RoyHub/nuclear/tests/examples/data/2xnr_mini.pdb'
crd_path = '/home/roy.gonzalez-aleman/RoyHub/nuclear/tests/examples/data/2xnr_RCXN010.crd'
mapping = {'H42': 100, "H5''": 95, 'H41': 72, "H1'": 65, "H2'": 50, "H2''": 26}

ligs = docking.CRD(crd_path).ag
rec = prd.parsePDB(pdb_path)
self = STDScorer(rec, ligs, mapping)

plt.scatter(range(len(self.scores)), self.scores)
plt.show()
