# Created by roy.gonzalez-aleman at 13/11/2023
"""
Control STDock workflow
"""
from os.path import basename, dirname, join

import matplotlib

matplotlib.use('agg')
import numpy as np
import numpy_indexed as npi
import pandas as pd
import prody as prd
import scipy.optimize as opt
from matplotlib import pyplot as plt
# from sklearn.metrics import r2_score
from tqdm import tqdm

import stdock.commons as cmn
import stdock.proj_paths as pp
from stdock.programs.lightdock import LightDock
from stdock.programs.vina import Vina
from stdock.rescorer import STDEpitope, STDScorer
from stdock.spectrum import Spectrum

prd.LOGGER.verbosity = 'none'


# todo: Move the generating plots section to the end and make it optional
# todo: Move plots and tables generation to another module
# todo: Uncomment the plotting section in the main function
# todo: Merge functions get_std0_data and get_std_af0_data
# todo: Check labels passed in std-epitopes / std-regions are unique (collectively, individual after splitting by - [1])

def std_func(x_data, std_max, k_sat):
    """
    Fit data to a mono-exponential curve

    Args:
        x_data: data in x-axis (saturation time)
        std_max: std max intensities
        k_sat: k saturation

    Returns:
        fitted data
    """
    return np.float64(std_max * (1 - np.exp(-k_sat * x_data)))


def kd_func(lig_conc, alpha, kd):
    """
    Fit data to calculate the dissociation constant

    Args:
        alpha: scaling factor
        lig_conc: ligand concentration
        kd: dissociation constant

    Returns:
        fitted data

    """
    return np.float64((alpha * lig_conc) / (kd + lig_conc))


def std_af_func(lig_conc, alpha, kd):
    """
    This function defines the equation for std_af(L) as a function of L, alpha,
     and Kd.

    Args:
      lig_conc: The independent variable, representing the ligand concentration.
      alpha: A scaling factor.
      kd: The dissociation constant.

    Returns:
      The calculated value of std_af(L) based on the given parameters.
    """
    return alpha * lig_conc / (lig_conc + kd)


def get_optimized_params(fit_func, x_data, y_data):
    """
    Get the optimized data
    Args:
        fit_func: function to which data will be fitted
        x_data: x-axis data (saturation times)
        y_data: y-axis data (STD intensities)

    Returns:
        info of fitting
    """
    return opt.curve_fit(fit_func, x_data, y_data, maxfev=50000)


class STDRunner:
    """
    Control STDock workflow
    """

    def __init__(self, config_args):
        """
        Args:
            config_args: arguments as parsed by STDConfig class
        """
        # Parse general arguments
        self.config_args = config_args

        # STD data
        self.std_data = None
        self.std_af_data = None
        self.kd = None
        self.kd_data = None

        # Spectra
        self.regions = None
        self.spectra = None
        self.intensities = None

        # Epitope mapping
        self.epitope = None
        self.mappings = None

        # Docking scores
        self.poses = None
        self.docking_obj = None
        self.docking_scores = None
        self.std_scores = None

    def run(self):
        """ Run STDock workflow step by step"""

        # [WF] MAP FROM SPECTRA
        if 'map-from-spectra' in self.config_args['template']:
            # Get intensities by integrating regions
            self.spectra = self.config_args['std-spectra']
            self.regions = self.config_args['std-regions']
            self.intensities = self.get_intensities()

            # Plot STD curves
            self.std_data = self.get_std0_data()
            self.plot_data()

            # Do epitope mapping (create files at STD directory)
            self.epitope = self.get_epitope_mappings(integrate=True)

            # Launch Kd determination
            lig_conc = self.epitope.keys()
            if (n_conc := len(lig_conc)) > 1:
                print(f'Multiple ({n_conc}) ligand concentrations detected.'
                      f' Launching Kd determination.')

                self.std_af_data = self.get_std_af0_data()
                self.kd_data, self.kd = self.get_kds()

        # [WF] MAP FROM VALUES
        if 'map-from-external' in self.config_args['template']:
            # Do epitope mapping
            self.epitope = self.config_args['std-epitopes']
            self.get_epitope_mappings()
            print('Epitope mapping completed')

        # [WF] THEN DOCK
        if 'dock' in self.config_args['template']:
            # Convert from pdb to pdbqt
            print('Launching docking as requested in config')
            self.config_args['ligand_pdbqt'] = cmn.pdb2pdbqt(
                self.config_args['ligand_pdb'], ligand_case=True)
            self.config_args['receptor_pdbqt'] = cmn.pdb2pdbqt(
                self.config_args['receptor_pdb'])

            # Docking of small molecules (THE ONLY CHOICE FOR PUBLIC RELEASE)
            if 'dock-small' in self.config_args['template']:
                self.docking_obj = self.launch_docking()
                topology, trajectory, rec_path = self.pre_stdscore_small()

                epitope = STDEpitope(self.config_args['STD'])
                qt_indices = self.reindex_pdb_qt()
                score_obj = STDScorer(epitope, topology=topology,
                                      trajectory=trajectory, rec_path=rec_path,
                                      qt_indices=qt_indices)

                self.std_scores = score_obj.scores
                self.docking_scores = self.parse_docking_scores()

            # Docking of peptides (ONLY FOR INTERNAL USE)
            elif 'dock-pept' in self.config_args['template']:
                # Prepare inputs for rescoring
                lig_str = f'lightdock_{basename(self.docking_obj.lig_path)}'
                rec_str = f'lightdock_{basename(self.docking_obj.rec_path)}'
                lig_path = next(cmn.recursive_finder(lig_str,
                                                     self.docking_obj.out_dir))
                rec_path = next(cmn.recursive_finder(rec_str,
                                                     self.docking_obj.out_dir))
                poses_dcd = self.docking_obj.poses_dcd

                # Compute std-score
                input_dir = self.config_args['output_dir']
                epitopes = STDEpitope(input_dir)
                score_obj = STDScorer(lig_path, rec_path, poses_dcd, epitopes)
                self.std_scores = score_obj.scores
                self.docking_scores = self.docking_obj.light_dock_scores

            else:
                raise ValueError('STDScore not implemented for this template')

            # Save scores
            self.save_scores()

    def get_intensities(self):
        """
        Get intensity values using TopSpin and the [std-regions] sections
        defined in the config file

        Returns:
            intensities: dict of intensities
        """
        spectra, regions = self.spectra, self.regions
        intensities = {}
        for lig_conc in tqdm(spectra, desc='Integrating Spectral Regions'):
            intensities.update({lig_conc: {}})

            for d20 in spectra[lig_conc]:
                off_proton = Spectrum(spectra[lig_conc][d20]['off'], regions)
                off_int = off_proton.integrate()
                off_int.columns = ['absolute']

                if 'on' in spectra[lig_conc][d20]:
                    on_proton = Spectrum(spectra[lig_conc][d20]['on'], regions)
                    on_int = on_proton.integrate()
                    on_int.columns = ['absolute']
                    std_intensities = (off_int - on_int) / off_int
                    intensities[lig_conc].update({d20: std_intensities})

                if 'diff' in spectra[lig_conc][d20]:
                    diff_proton = Spectrum(
                        spectra[lig_conc][d20]['diff'],
                        self.regions)
                    diff_int = diff_proton.integrate()
                    diff_int.columns = ['absolute']
                    std_intensities = diff_int / off_int
                    intensities[lig_conc].update({d20: std_intensities})
        return intensities

    def get_std_af0_data(self):
        """
        Reorder intensity data into per-H dataframes containing raw and fitted
        std-af data for each saturation time.

        Returns:
            reordered_data2: dict of raw and fitted dataframes
        """
        std_data = self.std_data
        prot_concentration = self.config_args['protein_conc']
        intensities = self.intensities
        first_key = list(intensities.keys())[0]
        t_sats = np.asarray(sorted(intensities[first_key]))
        use_exponential = False if len(t_sats) == 2 else True

        reordered_data2 = {}
        for lig_conc in std_data:
            reordered_data2.update({lig_conc: {}})
            for label in std_data[lig_conc]:

                # Create dataframe
                df = std_data[lig_conc][label]
                df['std_af'] = df['std'] * (lig_conc / prot_concentration)
                t_sats = df.tsat.to_numpy()
                stds_af = df.std_af.to_numpy()

                # Fit data to exponential curve
                if use_exponential:
                    try:
                        std_af_max, k_sat = \
                            get_optimized_params(std_func, t_sats,
                                                 stds_af)[0]
                        std_af_0 = std_af_max * k_sat
                    except RuntimeError:
                        print('Issues with fitting data. Setting values to 0')
                        std_af_max, k_sat, std_af_0 = 0, 0, 0
                    fit_data = std_func(t_sats, std_af_max, k_sat)
                    # Update dataframe
                    df['std_af_fit'] = fit_data
                    df['k_sat2'] = k_sat
                    df['std_af_max'] = std_af_max

                # Use the reduced-dataset approach
                else:
                    t_short = min(t_sats)
                    std_af_long = max(stds_af)
                    std_af_short = min(stds_af)
                    std_af_0 = (std_af_long / t_short) * np.log(
                        std_af_long / (std_af_long - std_af_short))
                df['std_af_0'] = std_af_0

                reordered_data2[lig_conc].update({label: df})
        if use_exponential:
            print('STD-AF data has been fitted to an exponential curve')
        else:
            print('STD-AF data derived through the reduced dataset approach')
        return reordered_data2

    def get_std0_data(self):
        """
        Reorder intensity data into per-H dataframes containing raw and fitted
        std data for each saturation time

        Returns:
            reordered_data: dict of raw and fitted dataframes
        """
        # Get hydrogen labels
        intensities = self.intensities
        first_key = list(intensities.keys())[0]
        second_key = list(intensities[first_key].keys())[0]
        labels = intensities[first_key][second_key].index.tolist()

        # Reorder data
        t_sats = np.asarray(sorted(intensities[first_key]))
        use_exponential = False if len(t_sats) == 2 else True

        reordered_data = {}
        for lig_conc in intensities:
            reordered_data.update({lig_conc: {}})
            for label in labels:

                # Create dataframe
                df = pd.DataFrame()
                df['tsat'] = t_sats
                stds = []
                for t_sat in t_sats:
                    stds.append(
                        intensities[lig_conc][t_sat].loc[label].absolute)
                df['std'] = stds

                # Fit data to exponential curve
                if use_exponential:
                    try:
                        std_max, k_sat = get_optimized_params(
                            std_func, t_sats, stds)[0]
                        std0 = std_max * k_sat
                    except RuntimeError:
                        print('Issues with fitting data. Setting values to 0')
                        std_max, k_sat, std0 = 0, 0, 0
                    fit_data = std_func(t_sats, std_max, k_sat)
                    df['std_fit'] = fit_data
                    df['k_sat'] = k_sat
                    df['std_max'] = std_max

                # Use the reduced-dataset approach
                else:
                    t_short = min(t_sats)
                    std_long = max(stds)
                    std_short = min(stds)
                    std0 = (std_long / t_short) * np.log(
                        std_long / (std_long - std_short))

                # Update dataframe
                df['std_0'] = std0
                reordered_data[lig_conc].update({label: df})

        if use_exponential:
            print('STD data has been fitted to an exponential curve')
        else:
            print('STD data derived through the reduced dataset approach')

        return reordered_data

    def plot_data(self):
        """
        Plot all STD buildup curves and fittings
        """
        # Styles
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        lines = ['-']
        shapes = ['o', 'v', 's', '*', 'd']
        styles = cmn.combinator([lines, shapes, colors])
        cmn.generic_matplotlib()

        # Plotting individual fitting data
        data = self.std_data
        for lig_conc in tqdm(data, desc=f'Plotting STD Data',
                             bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}"):
            for label in data[lig_conc]:
                df = data[lig_conc][label]
                self.plot_individual_data(f'{label}_{lig_conc}', df)

            # Plotting the ensemble of raw data
            self.plot_raw_ensemble(styles)

            # Plotting the ensemble of fitted data
            self.plot_fitted_ensemble(styles)

    def plot_individual_data(self, label, df):
        """
        Plot the raw data of each Hydrogen and its exponential fitting
        Args:
            label: a label for naming output plot
            df: dataframe containing data to plot
        """
        # Plot the actual data
        std_max = df.std_max[0]
        k_sat = df.k_sat[0]
        std0 = df.std_0[0]
        new_x = np.linspace(df.tsat.min(), df.tsat.max())

        plt.plot(df.tsat, df['std'], lw=0, label="Raw data", marker='x',
                 color='r',
                 alpha=0.5)
        plt.plot(df.tsat, df.std_fit, lw=0, label="Fitted data", marker='o',
                 color='navy')
        plt.plot(new_x, std_func(new_x, std_max, k_sat), lw=1, ls='--',
                 color='navy', alpha=0.5)
        plt.title(
            f'STD_max={std_max:.2e} | k_sat={k_sat:.2e} | STD0={std0:.2e}')
        plt.suptitle(f'{label}')
        plt.legend(loc='upper left')
        plt.xlabel('Saturation time (s)', fontweight='bold')
        plt.ylabel('STD intensity', fontweight='bold')
        out_name = join(self.config_args['STD'], f'fitted-{label}.png')
        plt.savefig(out_name, bbox_inches='tight', dpi=300)
        plt.close()

    def plot_raw_ensemble(self, styles):
        """
        Plot the raw data of all Hydrogens in the same graph
        Args:
            styles: matplotlib styles for legend construction
        """
        plt.xlabel('Saturation time (s)', fontweight='bold')
        plt.ylabel('STD intensity', fontweight='bold')

        data = self.std_data
        for lig_conc in data:
            for i, label in enumerate(data[lig_conc]):
                line, shape, color = styles[i]
                df = data[lig_conc][label]
                plt.plot(df.tsat, df['std'], lw=1, ls=line, marker=shape,
                         color=color,
                         ms=5, label=label, alpha=0.75)
            plt.legend(loc='right', bbox_to_anchor=(1.2, 0.5), ncol=1,
                       fontsize='small')
            out_name = join(self.config_args['STD'],
                            f'ENSEMBLE-RAW_{lig_conc}.png')
            plt.savefig(out_name, bbox_inches='tight', dpi=300)
            plt.close()

    def plot_fitted_ensemble(self, styles):
        """
        Plot the fitted data of all Hydrogens in the same graph
        Args:
            styles: matplotlib styles for legend construction
        """
        plt.xlabel('Saturation time (s)', fontweight='bold')
        plt.ylabel('STD intensity', fontweight='bold')

        data = self.std_data
        for lig_conc in data:
            for i, label in enumerate(data[lig_conc]):
                line, shape, color = styles[i]

                df = data[lig_conc][label]
                plt.plot(df.tsat, df.std_fit, lw=0, label=label, marker=shape,
                         color=color)

                std_max = df.std_max[0]
                k_sat = df.k_sat[0]
                new_x = np.linspace(df.tsat.min(), df.tsat.max())
                plt.plot(new_x, std_func(new_x, std_max, k_sat), lw=1, ls=line,
                         color=color, alpha=0.75)
            plt.legend(loc='right', bbox_to_anchor=(1.2, 0.5), ncol=1,
                       fontsize='small')
            out_name = join(self.config_args['STD'],
                            f'ENSEMBLE-FITTED_{lig_conc}.png')
            plt.savefig(out_name, bbox_inches='tight', dpi=300)
            plt.close()

    def get_epitope_mappings(self, integrate=False):
        """
        Get the epitope mapping for each ligand concentration

        Returns:
            mappings: a dict of {label: normalized std} for every ligand
                      concentration
        """
        # Get epitope mapping after spectra integration
        if integrate:
            data = self.std_data
            mappings = {}
            for lig_conc in data:
                mappings.update({lig_conc: {}})
                std0 = []
                labels = []
                max_val = 0
                for label in data[lig_conc]:
                    labels.append(label)
                    val = data[lig_conc][label].std_0[0]
                    std0.append(val)
                    if val > max_val:
                        max_val = val
                norm = np.asarray(std0) / max_val

                out_name = join(self.config_args['STD'],
                                f'EPITOPE-MAPPING_{lig_conc}.txt')
                with open(out_name, 'wt') as mapping:
                    for i, label in enumerate(labels):
                        mappings[lig_conc].update({label: norm[i]})
                        mapping.write(f'{label:>16}: {norm[i]:.2f}\n')

        # Get epitope mapping from external values specified in the config
        else:
            out_name = join(self.config_args['STD'],
                            f'EPITOPE-MAPPING_0.0.txt')
            labels = self.epitope.keys()
            norm = np.asarray(list(self.epitope.values()))
            mappings = {}
            with open(out_name, 'wt') as mapping:
                for i, label in enumerate(labels):
                    mappings.update({label: norm[i]})
                    mapping.write(f'{label:>16}: {norm[i]:.2f}\n')

        print('Epitope mapping completed')
        return mappings

    def reindex_pdb_qt(self):
        """
        Reindex pdb and pdbqt files to match the order of atoms.

        Returns:
            indices: indices that arranges pdb and pdbqt files atomic ordering
        """
        pdb_path = self.config_args['ligand_pdb']
        qt_path = self.config_args['ligand_pdbqt']
        pdb_parsed = prd.parsePDB(pdb_path)
        qt_parsed = cmn.Molecule(qt_path).parse()[0]
        indices = npi.indices(pdb_parsed.getCoords(), qt_parsed.getCoords(),
                              missing='mask')

        pdb_coords = pdb_parsed.getCoords()
        qt_coords = qt_parsed.getCoords()
        for coord in qt_coords:
            if coord not in pdb_coords:
                print('Missing coordinates in pdb file')

        return indices

    def launch_docking(self):
        """
        Launch docking engine job using [docking] section parameters in config
        """
        if 'dock-small' in self.config_args['template']:
            docking_obj = Vina(
                pp.vina_exe,
                self.config_args['receptor_pdbqt'],
                self.config_args['ligand_pdbqt'],
                self.config_args['num_poses'],
                self.config_args['rmsd_tolerance'],
                ['vina'],
                [self.config_args['exhaustiveness']],
                self.config_args['DOCKING'])
            docking_obj.run_commands()
            docking_obj.yield_filter_sort()

        elif 'dock-pept' in self.config_args['template']:
            docking_obj = LightDock(
                self.config_args['receptor_pdb'],
                self.config_args['ligand_pdb'],
                self.config_args['DOCKING'],
                num_steps=self.config_args['num_steps'],
                sf=self.config_args['scoring_function'])
        else:
            raise ValueError

        return docking_obj

    def save_scores(self):
        """
        Save std and docking scores
        """
        std_scores = self.std_scores
        docking_scores = self.docking_scores

        score_table = pd.DataFrame()
        sorted_keys = sorted(std_scores.keys())
        for lig_conc in sorted_keys:
            lig_scores = std_scores[lig_conc]
            if not lig_scores.size == len(docking_scores):
                raise ValueError(
                    'Mismatch in number of docking_scores & std scores')

            score_table[f'std_score_{lig_conc}'] = lig_scores

        score_table['docking_score'] = docking_scores
        out_name = join(self.docking_obj.out_dir, 'scores_all.txt')
        with open(out_name, 'wt') as scores:
            score_table.to_string(scores, index=False)

    def pre_stdscore_small(self):
        """
        Prepare inputs for getting scores in the case of a small docking run.

        Returns:
            topo_path: path to the topology file
            traj_path: path to the trajectory file
            rec_path: path to the receptor file
            epitope_path: path to the epitope mapping file
        """
        # find poses.pdb
        docking_dir = self.config_args['DOCKING']
        std_dir = self.config_args['STD']
        poses_path = next(cmn.recursive_finder('poses.pdb', docking_dir))
        poses_parsed = prd.parsePDB(poses_path)

        # Get dcd trajectory
        traj_path = join(dirname(docking_dir), 'trajectory.dcd')
        prd.writeDCD(traj_path, poses_parsed)

        # Get the pdb topology
        topo_path = join(dirname(docking_dir), 'topology.pdb')
        prd.writePDB(topo_path, poses_parsed, csets=0)

        # Get the receptor_path
        rec_pdbqt = self.config_args['receptor_pdbqt']
        rec_parsed = cmn.Molecule(rec_pdbqt).parse()[0]
        rec_path = join(dirname(docking_dir), 'receptor.pdb')
        prd.writePDB(rec_path, rec_parsed)

        # Get the epitope mapping path
        # epitope_path = join(std_dir, 'EPITOPE-MAPPING_0.5.txt')
        # with open(epitope_path, 'wt') as ep:
        #     for key, value in self.epitope.items():
        #         ep.write(f'{key} : {value}\n')

        return topo_path, traj_path, rec_path

    def parse_docking_scores(self):
        """
        Parse docking scores from the docking output directory

        Returns:
            scores: numpy array of docking scores
        """
        dock_out_dir = self.docking_obj.out_dir
        scores_path = next(
            cmn.recursive_finder('scores.txt', dock_out_dir))
        parsed = pd.read_table(scores_path, sep='\s+', header=None).T.iloc[1]
        scores = np.asarray(parsed)
        return scores

    def get_kds(self):
        """
        Fit the data to get the dissociation constant.

        Returns:
            reordered_data: dict of reordered data
            kd: dissociation constant
        """
        data = self.std_af_data
        lig_concs = list(data.keys())
        labels = list(data[lig_concs[0]].keys())

        reordered_data = {}
        for label in labels:
            reordered_data.update({label: {'conc': [], 'std_af0': []}})
            for lig_conc in lig_concs:
                reordered_data[label]['conc'].append(lig_conc)
                std_af0 = data[lig_conc][label].std_af_0[0]
                reordered_data[label]['std_af0'].append(std_af0)

            # Fit data
            try:
                alpha, kd = get_optimized_params(kd_func,
                                                 reordered_data[label]['conc'],
                                                 reordered_data[label][
                                                     'std_af0'])[0]
            except RuntimeError:
                print('Issues with fitting data. Setting values to 0')
                alpha, kd = 0, 0

            # Update dataframe
            reordered_data[label]['alpha'] = alpha
            reordered_data[label]['kd'] = kd

        # Get kd as the average of all positive kd values
        values = np.asarray(
            [reordered_data[label]['kd'] for label in reordered_data])
        positive = values[values > 0]
        kd = np.abs(positive).mean()
        print('Kd determination completed')
        return reordered_data, kd


# =============================================================================
# Debugging Case studies
# =============================================================================

# config_path = '/home/gonzalezroy/RoyHub/stdock/tests/paper/hur/map-from-values-then-dock_hur.cfg'
# config_path = '/home/gonzalezroy/RoyHub/stdock/tests/paper/imaging/config_kd.cfg'
# config_path = '/home/gonzalezroy/RoyHub/stdock/tests/paper/reduced_dataset/config_rd.cfg'
# config_path = '/home/gonzalezroy/RoyHub/stdock/tests/paper/imaging-docking/config_kd_docking.cfg'


# config_path = '/home/gonzalezroy/RoyHub/stdock/tests/paper/hur/HuR_reproduction/M2_e/map-from-values-then-dock_hur.cfg'
# config_path = '/home/gonzalezroy/RoyHub/stdock/tests/paper/hur/HuR_reproduction/M3_e/map-from-values-then-dock_hur.cfg'
# config_path = '/home/gonzalezroy/RoyHub/stdock/tests/paper/hur/HuR_reproduction/M4_e/map-from-values-then-dock_hur.cfg'
# config_path = '/home/gonzalezroy/RoyHub/stdock/tests/paper/hur/HuR_reproduction/M6_e/map-from-values-then-dock_hur.cfg'
# config_path = '/home/gonzalezroy/RoyHub/stdock/tests/paper/hur/HuR_reproduction/M7_e/map-from-values-then-dock_hur.cfg'
# config_path = '/home/gonzalezroy/RoyHub/stdock/tests/paper/hur/HuR_reproduction/M11_e1/map-from-values-then-dock_hur.cfg'
# config_path = '/home/gonzalezroy/RoyHub/stdock/tests/paper/hur/HuR_reproduction/M12_e/map-from-values-then-dock_hur.cfg'
# config_path = '/home/gonzalezroy/RoyHub/stdock/tests/paper/error/map-from-spectra-then-dock-small.cfg'

configs = {
    'map-from-spectra': '/home/gonzalezroy/RoyHub/stdock/tests/paper/trolls/map-from-spectra.cfg',
    'map-from-spectra-then-dock-small': '/home/gonzalezroy/RoyHub/stdock/tests/paper/error/map-from-spectra-then-dock-small.cfg',
    'map-from-values-then-dock': '/home/gonzalezroy/RoyHub/stdock/tests/paper/hur/HuR_reproduction/M1/map-from-values-then-dock_hur.cfg',
    # 'map-from-spectra-then-dock-pept': 'map-from-spectra-then-dock-pept.cfg',
}

# config_path = configs['map-from-spectra-then-dock-small']
# config_path = configs['map-from-spectra']
# config_path = configs['map-from-values-then-dock']
# params = cfg.allowed_parameters
# valid_templates = cfg.allowed_templates
# args = cfg.STDConfig(config_path, params, valid_templates).config_args
# self = STDRunner(args)
# self.run()

# =============================================================================
#
# =============================================================================
# import numpy as np
#
#
# # Example usage:
# L_data = np.array(
#     [0.2, 0.4, 0.6, 0.8, 1]) * 1000  # Example ligand concentrations
# std_af_data = np.array([4.6, 5.5, 7.0, 7.1, 7.2])  # Example std_af(L) data
#
# popt, pcov = get_optimized_params(kd_func, L_data, std_af_data)
#
# # Extract the optimized values
# alpha_opt = popt[0]
# Kd_opt = popt[1]
#
# print(f"Optimized alpha: {alpha_opt}")
# print(f"Optimized Kd: {Kd_opt}")
#
# # Extract the optimized values
# alpha_opt = popt[0]
# Kd_opt = popt[1]
#
# # Calculate the standard errors
# std_err_alpha = np.sqrt(pcov[0, 0])
# std_err_Kd = np.sqrt(pcov[1, 1])
#
#
# # Round the values to the nearest integer (no decimal places)
# Kd_rounded = round(Kd_opt)
# std_err_Kd_rounded = round(std_err_Kd)
#
# # Report the result
# print(f"Kd: {Kd_rounded} ± {std_err_Kd_rounded}")
