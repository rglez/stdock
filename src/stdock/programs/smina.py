# Created by roy.gonzalez-aleman at 03/01/2024
import re
from os.path import join, split

import prody as prd
from scipy.spatial import cKDTree as ckd

import commons as cmn
from programs.vina import Vina


class Smina(Vina):
    def get_commands(self):
        commands = []
        for exhaustiveness in self.exhaustiveness_list:
            for sf in self.scoring_functions:
                sub_dir = f'smina_{sf}_{cmn.syn[exhaustiveness]}'
                out_dir = join(self.out_dir, 'smina', sub_dir)
                cmd = (f'{self.exe}'
                       f' --receptor {self.rec_path}'
                       f' --ligand {self.lig_path}'
                       f' --scoring {sf}'
                       f' --autobox_ligand {self.rec_path}'
                       f' --seed 910623'
                       f' --exhaustiveness {exhaustiveness}'
                       f' --num_modes {self.num_poses}'
                       f' --min_rmsd_filter {self.rmsd_tol}'
                       f' --energy_range 10000'
                       f' --out {join(out_dir, sub_dir)}.pdbqt'
                       )
                commands.append(cmd)
        return commands, None

    def yield_filter_sort(self):
        string_score = 'REMARK minimizedAffinity.*'

        for pdbqt in cmn.recursive_finder('*pdbqt', self.out_dir):
            outdir = split(pdbqt)[0]

            # Get all scores
            with open(pdbqt, 'rt') as inp:
                file_string = inp.read()
            score_raw = re.findall(string_score, file_string)
            scores = [x.split()[-1].strip() for x in score_raw]
            scores_path = join(outdir, 'scores.txt')
            with open(scores_path, 'wt') as out:
                for i, x in enumerate(scores):
                    out.write(f'{i}    {x}\n')

            # Get all poses
            ensemble = cmn.Molecule(pdbqt).get_ensemble()
            out_name = join(outdir, 'poses.pdb')
            prd.writePDB(out_name, ensemble)

            # Get filtered indices
            rec_kdt = ckd(self.rec_parsed.getCoords())
            lig_parsed = cmn.Molecule(pdbqt).parse()
            filtered_indices, filtered_ligs = cmn.get_filtered_indices(
                rec_kdt, lig_parsed)

            # Get filtered poses
            lig_filtered = prd.Ensemble()
            [lig_filtered.addCoordset(x.getCoords()) for x in filtered_ligs]
            lig_filtered.setAtoms(filtered_ligs[0])
            prd.writePDB(join(outdir, 'poses_filtered.pdb'), lig_filtered)

            # Get filtered scores
            with open(join(outdir, 'scores_filtered.txt'), 'wt') as out:
                for index in filtered_indices:
                    out.write(f'{index}    {scores[index]}\n')
