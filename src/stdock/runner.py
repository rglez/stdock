# Created by roy.gonzalez-aleman at 13/11/2023
from os.path import join, basename, dirname

import numpy as np
import numpy_indexed as npi
import pandas as pd
import prody as prd
import scipy.optimize as opt
from bruker.api.topspin import Topspin
from bruker.data.nmr import *
from matplotlib import pyplot as plt

import commons as cmn
import proj_paths as pp
from programs.lightdock import LightDock
from programs.vina import Vina
from rescorer import STDEpitope, STDScorer

prd.LOGGER.verbosity = 'none'


def fit_func(x_data, std_max, k_sat):
    """
    Fit data to a mono-exponential curve

    Args:
        x_data: data in x-axis (saturation time)
        std_max: std max intensities
        k_sat: k saturation

    Returns:
        fitted data
    """
    return std_max * (1 - np.exp(-k_sat * x_data))


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
    fit_info = opt.curve_fit(fit_func, x_data, y_data, maxfev=50000)
    return fit_info[0]


def get_integral(proton, start, end):
    """
    Get the absolute integral value of a spectrum region

    Args:
        proton: proton object from TopSpin
        start: start value for integration
        end: end value for integration

    Returns:
        summation: absolute integral value
    """
    res = proton.getSpecDataPoints()
    if res.get(EXCEPTION):
        return -1

    startIndex = proton.getIndexFromPhysical(start, 0)
    endIndex = proton.getIndexFromPhysical(end, 0)

    summation = 0
    data = res.get(DATA_POINTS)
    for i in range(startIndex, endIndex):
        summation = summation + data[i]
    return summation


class Spectrum:
    """
    Deal with spectrum-related tasks
    """

    def __init__(self, pdata_path, regions, out_dir='.'):
        # Pdata path & sub folders
        self.path = cmn.check_path(pdata_path)
        self.regions = regions

        # Out dir
        self.root_id = os.sep.join(self.path.split(os.sep)[-3:-1])
        self.out_dir = join(out_dir, self.root_id)

        # TopSpin interface
        self.top = Topspin()
        self.dp = self.top.getDataProvider()

    def integrate(self):
        """
        Integrates the spectrum in the defined regions

        Returns:
            absolute: dataframe of absolute integrals
        """
        # Get TopSpin proton object
        integrals = cmn.recursive_defaultdict()
        proton = self.dp.getNMRData(self.path)
        if proton is None:
            raise ValueError(
                'Something went wrong :( Check the Top Spin connection.')

        # Integrate each region defined in the config file
        for index in self.regions:
            start = self.regions[index]['upper']
            end = self.regions[index]['lower']

            absolute_integral = get_integral(proton, start, end)
            integrals[self.path][index] = absolute_integral

        # Formatting absolute & relative integrals
        absolute = pd.DataFrame(integrals).astype(int)
        absolute['labels'] = [self.regions[x]['label'] for x in self.regions]
        absolute.set_index(absolute.labels, inplace=True)
        del (absolute['labels'])
        return absolute


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

        # Automatic conversion of pdb input files to pdbqt
        self.config_args['ligand_pdbqt'] = cmn.pdb2pdbqt(
            config_args['ligand_pdb'], ligand_case=True)

        self.config_args['receptor_pdbqt'] = cmn.pdb2pdbqt(
            config_args['receptor_pdb'])

        # Conditional parsing of sections
        if 'map-from-spectra' in self.config_args['template']:
            self.spectra = self.config_args['std-spectra']
            self.regions = self.config_args['std-regions']

        if 'map-from-values' in self.config_args['template']:
            self.epitope = self.config_args['std-epitopes']

        # To be filled during running
        self.intensities = None
        self.data = None
        self.mappings = None
        self.poses = None
        self.std_scores = None
        self.docking_obj = None
        self.docking_scores = None

    def run(self):
        """
        Run STDock workflow step by step
        """

        # Map from spectra
        if 'map-from-spectra' in self.config_args['template']:
            # 1. Get intensities by integrating regions
            print('Starting integrations')
            self.intensities = self.get_intensities()
            print('Declared regions have been integrated\n')

            # 2. Plot STD curves
            self.data = self.reorder_data()
            print('Starting to plot')
            self.plot_data()
            print('Plots have been produced')

            # 3. Do epitope mapping
            self.mappings = self.get_epitope_mappings()
            print('Epitope mapping completed')

        # Map from external values
        if 'map-from-values' in self.config_args['template']:
            # 3. Do epitope mapping
            self.mappings = self.epitope
            print('Epitope mapping completed')

        # Docking
        if 'dock' in self.config_args['template']:
            # 4. Launch docking
            print('Launching docking as requested in config')
            self.docking_obj = self.launch_docking()

            # 5. Rescore using stdscore
            input_dir = self.config_args['output_dir']

            if 'dock-pept' in self.config_args['template']:
                # 5.1 Prepare inputs for rescoring
                lig_str = f'lightdock_{basename(self.docking_obj.lig_path)}'
                rec_str = f'lightdock_{basename(self.docking_obj.rec_path)}'
                lig_path = next(cmn.recursive_finder(lig_str,
                                                     self.docking_obj.out_dir))
                rec_path = next(cmn.recursive_finder(rec_str,
                                                     self.docking_obj.out_dir))
                poses_dcd = self.docking_obj.poses_dcd

                # 5.2 Compute std-score
                epitopes = STDEpitope(input_dir)
                score_obj = STDScorer(lig_path, rec_path, poses_dcd, epitopes)
                self.std_scores = score_obj.scores
                self.docking_scores = self.docking_obj.light_dock_scores

            elif 'dock-small' in self.config_args['template']:
                topology, trajectory, rec_path, epitope_raw = self.pre_stdscore_small()
                qt_indices = self.reindex_pdb_qt()
                epitope = STDEpitope(dirname(epitope_raw))
                score_obj = STDScorer(epitope, topology=topology,
                                      trajectory=trajectory, rec_path=rec_path,
                                      qt_indices=qt_indices)
                self.std_scores = score_obj.scores
                self.docking_scores = self.parse_docking_scores()
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
        for lig_conc in spectra:
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

    def reorder_data(self):
        """
        Reorder intensity data into per-H dataframes containing raw and fitted
        data for each saturation time

        Returns:
            reordered_data: dict of raw and fitted dataframes

        """
        intensities = self.intensities

        # Hydrogen labels
        first_key = list(intensities.keys())[0]
        second_key = list(intensities[first_key].keys())[0]
        labels = intensities[first_key][second_key].index.tolist()

        # Reorder data
        columns = ['x', 'y', 'y_fit', 'k_sat', 'std_max', 'std_0']
        x_data = np.asarray(sorted(intensities[first_key]))

        reordered_data = {}
        for lig_conc in intensities:
            reordered_data.update({lig_conc: {}})
            for label in labels:
                df = pd.DataFrame(columns=columns)
                df.x = x_data
                y_data = []
                for d20 in x_data:
                    y_data.append(
                        intensities[lig_conc][d20].loc[label].absolute)
                df.y = y_data
                try:
                    std_max, k_sat = get_optimized_params(
                        fit_func, x_data, y_data)
                    std0 = std_max * k_sat
                except RuntimeError:
                    std_max, k_sat, std0 = 0, 0, 0
                fit_data = fit_func(x_data, std_max, k_sat)
                df.y_fit = fit_data
                df.k_sat = k_sat
                df.std_max = std_max
                df.std_0 = std0
                reordered_data[lig_conc].update({label: df})
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
        data = self.data
        for lig_conc in data:
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
        new_x = np.linspace(df.x.min(), df.x.max())

        plt.plot(df.x, df.y, lw=0, label="Raw data", marker='x', color='r',
                 alpha=0.5)
        plt.plot(df.x, df.y_fit, lw=0, label="Fitted data", marker='o',
                 color='navy')
        plt.plot(new_x, fit_func(new_x, std_max, k_sat), lw=1, ls='--',
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

        data = self.data
        for lig_conc in data:
            for i, label in enumerate(data[lig_conc]):
                line, shape, color = styles[i]
                df = data[lig_conc][label]
                plt.plot(df.x, df.y, lw=1, ls=line, marker=shape, color=color,
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

        data = self.data
        for lig_conc in data:
            for i, label in enumerate(data[lig_conc]):
                line, shape, color = styles[i]

                df = data[lig_conc][label]
                plt.plot(df.x, df.y_fit, lw=0, label=label, marker=shape,
                         color=color)

                std_max = df.std_max[0]
                k_sat = df.k_sat[0]
                new_x = np.linspace(df.x.min(), df.x.max())
                plt.plot(new_x, fit_func(new_x, std_max, k_sat), lw=1, ls=line,
                         color=color, alpha=0.75)
            plt.legend(loc='right', bbox_to_anchor=(1.2, 0.5), ncol=1,
                       fontsize='small')
            out_name = join(self.config_args['STD'],
                            f'ENSEMBLE-FITTED_{lig_conc}.png')
            plt.savefig(out_name, bbox_inches='tight', dpi=300)
            plt.close()

    def get_epitope_mappings(self):
        """
        Get the epitope mapping for each ligand concentration

        Returns:
            mappings: a dict of {label: normalized std} for every ligand
                      concentration
        """
        data = self.data
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
        return mappings

    def reindex_pdb_qt(self):
        pdb_path = self.config_args['ligand_pdb']
        qt_path = self.config_args['ligand_pdbqt']
        pdb_parsed = prd.parsePDB(pdb_path)
        qt_parsed = cmn.Molecule(qt_path).parse()[0]
        indices = npi.indices(pdb_parsed.getCoords(), qt_parsed.getCoords())
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
        for lig_conc in std_scores:
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
        epitope_path = join(std_dir, 'EPITOPE-MAPPING_0.5.txt')
        with open(epitope_path, 'wt') as ep:
            for key, value in self.epitope.items():
                ep.write(f'{key} : {value}\n')

        return topo_path, traj_path, rec_path, epitope_path

    def parse_docking_scores(self):
        dock_out_dir = self.docking_obj.out_dir
        scores_path = next(cmn.recursive_finder('scores_filtered.txt', dock_out_dir))
        parsed = pd.read_table(scores_path, sep='\s+', header=None).T.iloc[1]
        return np.asarray(parsed)


# =============================================================================
# Debugging data
# =============================================================================
# #### Debugging Spectrum
# pdata_path = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/STDCK2a-300v2/12/pdata/11/'
# regions = config_args['std-regions']
# self = Spectrum(pdata_path, regions)

# , total=len(self.commands)
# =============================================================================
# Debugging HUr
# =============================================================================
# import config as cfg
#
# config_path = '/home/gonzalezroy/RoyHub/stdock/tests/paper/hur/map-from-values-then-dock_hur.cfg'
# params = cfg.allowed_parameters
# valid_templates = cfg.allowed_templates
# args = cfg.STDConfig(config_path, params, valid_templates).config_args
# self = STDRunner(args)
# self.run()
#
# import matplotlib.pyplot as plt
# scores = pd.read_table('/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/00-CASE_STUDY/HuR/M9/stdock-M9/DOCKING/vina_vina_medium/scores.txt',
#                        sep='\s+', header=None)
# plt.scatter(self.std_scores[0.5], scores[1])
