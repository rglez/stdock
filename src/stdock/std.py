# Created by roy.gonzalez-aleman at 13/11/2023
from os.path import join

import matplotlib as mpl
import numpy as np
import pandas as pd
import scipy.optimize as opt
from bruker.api.topspin import Topspin
from bruker.data.nmr import *
from matplotlib import pyplot as plt

import commons as cmn


def fit_func(x_data, std_max, k_sat):
    return std_max * (1 - np.exp(-k_sat * x_data))


def get_optimized_params(fit_func, x_data, y_data):
    fit_info = opt.curve_fit(fit_func, x_data, y_data, maxfev=50000)
    return fit_info[0]


def generic_matplotlib():
    mpl.rc('figure', figsize=[12, 8], dpi=300)
    mpl.rc('xtick', direction='in', top=True)
    mpl.rc('xtick.major', top=False, )
    mpl.rc('xtick.minor', top=True, visible=True)
    mpl.rc('ytick', direction='in', right=True)
    mpl.rc('ytick.major', right=True, )
    mpl.rc('ytick.minor', right=True, visible=True)

    mpl.rc('axes', labelsize=24)
    mpl.rc('lines', linewidth=12, color='k')
    mpl.rc('font', family='monospace', size=16)
    mpl.rc('grid', alpha=0.5, color='gray', linewidth=1, linestyle='--')


def reset_matplotlib():
    mpl.rcParams.update(mpl.rcParamsDefault)


def get_integral(proton, start, end):
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
    def __init__(self, pdata_path, regions, out_dir='.'):
        # Pdata path & sub folders
        self.path = pdata_path
        cmn.check_path(self.path)
        self.regions = regions

        # Out dir
        self.root_id = os.sep.join(self.path.split(os.sep)[-3:-1])
        self.out_dir = join(out_dir, self.root_id)

        # TopSpin interface
        self.top = Topspin()
        self.dp = self.top.getDataProvider()

    def integrate(self):
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

    def save_integrations(self, integrations_list, labels):
        for i, integrations in enumerate(integrations_list):
            if integrations is not None:
                log_path = join(self.out_dir, f'integrations_{labels[i]}.txt')
                os.makedirs(self.out_dir, exist_ok=True)
                with open(log_path, 'wt') as log:
                    integrations.to_string(log)


class STDRunner:
    def __init__(self, config_args):
        self.config_args = config_args
        self.spectra = self.config_args['std-spectra']

        self.regions = self.config_args['std-regions']
        print('Starting integrations')

        self.intensities = self.get_intensities()
        print('Declared regions have been integrated')

        self.data = self.reorder_data()
        print('Starting to plot')

        self.plot_data()
        print('Plots have been produced')

        self.do_epitope_mapping()
        print('Epitope mapping completed')

    # todo: refactor this function
    def get_intensities(self):
        job_type = self.config_args['integration_type']
        intensities = {}
        if job_type == 2:
            for d20 in self.spectra:
                off_proton = Spectrum(self.spectra[d20]['off'], self.regions)
                off_integrals = off_proton.integrate()
                off_integrals.columns = ['absolute']

                if 'on' in self.spectra[d20]:
                    on_proton = Spectrum(self.spectra[d20]['on'], self.regions)
                    on_integrals = on_proton.integrate()
                    on_integrals.columns = ['absolute']
                    numerator = off_integrals - on_integrals
                    std_intensities = numerator / off_integrals
                    intensities.update({d20: std_intensities})

                if 'diff' in self.spectra[d20]:
                    diff_proton = Spectrum(self.spectra[d20]['diff'],
                                           self.regions)
                    diff_integrals = diff_proton.integrate()
                    diff_integrals.columns = ['absolute']
                    std_intensities = diff_integrals / off_integrals
                    intensities.update({d20: std_intensities})

        elif job_type == 3:
            for lig_conc in self.spectra:
                self.spectra.update({lig_conc: {}})

                for d20 in self.spectra[lig_conc]:
                    off_proton = Spectrum(self.spectra[d20]['off'],
                                          self.regions)
                    off_integrals = off_proton.integrate()
                    off_integrals.columns = ['absolute']

                    if 'on' in self.spectra[d20]:
                        on_proton = Spectrum(self.spectra[d20]['on'],
                                             self.regions)
                        on_integrals = on_proton.integrate()
                        on_integrals.columns = ['absolute']
                        std_intensities = (
                                                  off_integrals - on_integrals) / off_integrals
                        intensities[lig_conc].update({d20: std_intensities})

                    if 'off' in self.spectra[d20]:
                        diff_proton = Spectrum(self.spectra[d20]['diff'],
                                               self.regions)
                        diff_integrals = diff_proton.integrate()
                        diff_integrals.columns = ['absolute']
                        std_intensities = diff_integrals / off_integrals
                        intensities[lig_conc].update({d20: std_intensities})
        else:
            raise ValueError('The job type must be  2 or 3')
        return intensities

    def reorder_data(self):
        # Hydrogen labels
        labels = list(self.intensities.values())[0].index.tolist()
        columns = ['x', 'y', 'y_fit', 'k_sat', 'std_max', 'std_0']
        x_data = np.asarray(sorted(self.intensities.keys()))

        # Reorder data
        reordered_data = {}
        for label in labels:
            df = pd.DataFrame(columns=columns)
            df.x = x_data
            y_data = []
            for d20 in x_data:
                y_data.append(self.intensities[d20].loc[label].absolute)
            df.y = y_data
            try:
                std_max, k_sat = get_optimized_params(fit_func, x_data, y_data)
                std0 = std_max * k_sat
            except RuntimeError:
                std_max, k_sat, std0 = 0, 0, 0
            fit_data = fit_func(x_data, std_max, k_sat)
            df.y_fit = fit_data
            df.k_sat = k_sat
            df.std_max = std_max
            df.std_0 = std0
            reordered_data.update({label: df})
        return reordered_data

    def plot_data(self):
        # Styles
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        lines = ['-']
        shapes = ['o', 'v', 's', '*', 'd']
        styles = cmn.combinator([lines, shapes, colors])
        generic_matplotlib()

        # Plotting individual fitting data
        for label in self.data:
            df = self.data[label]
            self.plot_individual_data(label, df)

        # Plotting the ensemble of raw data
        self.plot_raw_ensemble(styles)

        # Plotting the ensemble of fitted data
        self.plot_fitted_ensemble(styles)

    def plot_individual_data(self, label, df):
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
        plt.xlabel('Saturation time (s)', fontweight='bold')
        plt.ylabel('STD intensity', fontweight='bold')
        for i, label in enumerate(self.data):
            line, shape, color = styles[i]
            df = self.data[label]
            plt.plot(df.x, df.y, lw=1, ls=line, marker=shape, color=color,
                     ms=5, label=label, alpha=0.75)
        plt.legend(loc='right', bbox_to_anchor=(1.2, 0.5), ncol=1,
                   fontsize='small')
        out_name = join(self.config_args['STD'], 'ENSEMBLE-RAW.png')
        plt.savefig(out_name, bbox_inches='tight', dpi=300)
        plt.close()

    def plot_fitted_ensemble(self, styles):
        plt.xlabel('Saturation time (s)', fontweight='bold')
        plt.ylabel('STD intensity', fontweight='bold')
        for i, label in enumerate(self.data):
            line, shape, color = styles[i]
            df = self.data[label]
            std_max = df.std_max[0]
            k_sat = df.k_sat[0]
            new_x = np.linspace(df.x.min(), df.x.max())
            plt.plot(df.x, df.y_fit, lw=0, label=label, marker=shape,
                     color=color)
            plt.plot(new_x, fit_func(new_x, std_max, k_sat), lw=1, ls=line,
                     color=color, alpha=0.75)
        plt.legend(loc='right', bbox_to_anchor=(1.2, 0.5), ncol=1,
                   fontsize='small')
        out_name = join(self.config_args['STD'], 'ENSEMBLE-FITTED.png')
        plt.savefig(out_name, bbox_inches='tight', dpi=300)
        plt.close()

    def do_epitope_mapping(self):
        std0 = []
        labels = []
        max_val = 0
        for label in self.data:
            labels.append(label)
            val = self.data[label].std_0[0]
            std0.append(val)
            if val > max_val:
                max_val = val
        norm = np.asarray(std0) / max_val
        out_name = join(self.config_args['STD'], 'EPITOPE-MAPPING.txt')
        with open(out_name, 'wt') as mapping:
            for i, label in enumerate(labels):
                mapping.write(f'{label:>16}: {norm[i]:.2f}\n')


# =============================================================================
# Debugging data
# =============================================================================

# #### Debugging Spectrum
# pdata_path = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/STDCK2a-300v2/12/pdata/11/'
# regions = config_args['std-regions']
# self = Spectrum(pdata_path, regions)


# #### Debugging STDRunner
# import config as cfg
#
# config_path = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/config.cfg'
# params = cfg.allowed_parameters
# valid_templates = cfg.allowed_templates
# args = cfg.STDConfig(config_path, params, valid_templates).config_args
# self = STDRunner(args)
