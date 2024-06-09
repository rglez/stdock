# Created by roy.gonzalez-aleman at 26/02/2024
import commons as cmn
import protocols


class Workflow:
    """
    Run workflow to elucidate CIGB300-CK2_ALPHA binding modes
    """

    def __init__(self, rec_init_path, lig_init_path):
        self.rec_init = cmn.check_path(rec_init_path)
        self.lig_init = cmn.check_path(lig_init_path)

    def run(self):
        """
        Runs the workflow
        """
        # 1. MINIMIZE RECEPTOR
        rec_minim = protocols.minimize_receptor(self.rec_init)
        # 2. GENERATE LIG CONFORMERS
        # 2.1 CLUSTERING LIG CONFORMERS
        # 2.2 MINIMIZE CONFORMERS SEED
        # 3. RUN DOCKING FOR EACH CONFORMER

        # 4. MINIMIZE REC-LIG COMPLEXES

        # COMPUTE STDSCORE
