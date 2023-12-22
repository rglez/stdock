# Created by roy.gonzalez-aleman at 13/11/2023
import configparser
import os
from os.path import abspath, dirname, isabs, join, normpath

import numpy as np

import commons as cmn

#: Total number of cores
num_cores = os.cpu_count()

#: Allowed section templates in the config file
allowed_templates = {
    'map-from-spectra':
        {'generals', 'std-regions', 'std-spectra'},

    'map-from-values-then-dock':
        {},

    'map-from-spectra-then-dock':
        {},

    'dock-without-map':
        {}
}

#: Allowed keys in the config file (dtypes & expected values)
allowed_parameters = {

    # General parameters
    'generals': {
        'output_dir': {'dtype': 'path', 'check_exist': False}},

    # STD-NMR-related parameters
    'std-regions': None,

    # (On-res, off-res, spectra) tuples
    'std-spectra': None
}


class Param:
    def __init__(self, key, value, *args, **kwargs):
        self.key = key
        self.value = value
        self.args = args
        self.kwargs = kwargs

    def check(self):
        raise NotImplementedError()


class NumericParam(Param):
    def check(self):
        dtype = self.kwargs['dtype']
        minim = self.kwargs['min']
        maxim = self.kwargs['max']
        cmn.check_numeric_in_range(self.key, self.value, dtype, minim, maxim)


class PathParam(Param):
    def check(self):
        path = self.value
        cmn.check_path(path, check_exist=self.kwargs['check_exist'])


class ChoiceParam(Param):
    def check(self):
        choices = self.kwargs['values']
        if choices and (not (self.value in choices)):
            raise ValueError(
                f'\n Error in {self.key}. Passed "{self.value}" but available'
                f' options are: {choices}.')


class Config:
    def __init__(self, config_path, legal_params, legal_templates):

        # Parse class args
        self.config_path = cmn.check_path(config_path, check_exist=True)
        self.legal_params = legal_params
        self.legal_templates = legal_templates

        # Parsing from class args
        self.config_dir = abspath(dirname(self.config_path))
        self.keyless_sections = self.detect_keyless_sections()
        self.config_obj = self.read_config_file()
        self.template = self.detect_template()

        # Run checks
        self.check_missing_keys()
        self.config_args = self.check_params()
        self.check_constraints()
        self.config_args['template'] = self.template

    def detect_keyless_sections(self):
        params = self.legal_params
        return [x for x in params if params[x] is None]

    def read_config_file(self):
        config_obj = configparser.ConfigParser(allow_no_value=True,
                                               inline_comment_prefixes='#')
        config_obj.optionxform = str
        config_obj.read(self.config_path)
        return config_obj

    def detect_template(self):
        current_config = set(self.config_obj.sections())
        for template_name, template in self.legal_templates.items():
            if template == current_config:
                return template_name
        raise ValueError(
            f'Declared sections in the configuration file do not '
            f'correspond to any pre-defined template. Currently allowed '
            f'templates are: {self.legal_templates}.')

    def check_missing_keys(self):
        current_template = self.legal_templates[self.template].copy()
        current_params = self.legal_params
        [current_template.remove(x) for x in self.keyless_sections]

        for section in current_template:
            config_file_keys = list(self.config_obj[section].keys())
            for key in current_params[section]:
                if key not in config_file_keys:
                    raise KeyError(
                        f'Key "{key}" is missing from the config'
                        f' file. Please specify its value.')

    def check_params(self):
        config_args = dict()

        config_dir = self.config_dir
        parsed_sections = self.config_obj.sections().copy()
        [parsed_sections.remove(x) for x in self.keyless_sections]

        for section in parsed_sections:
            items = self.config_obj[section].items()
            for key, value in items:

                param_info = self.legal_params[section][key]
                dtype = param_info['dtype']
                if dtype in {float, int}:
                    param_obj = NumericParam(key, dtype(value), **param_info)
                elif dtype == 'path':
                    value = value if isabs(value) else join(config_dir, value)
                    param_obj = PathParam(key, value, **param_info)
                elif dtype == str:
                    param_obj = ChoiceParam(key, value, **param_info)
                else:
                    raise ValueError(
                        f"\n{section}.{key}'s dtype is wrong: {dtype}.")
                param_obj.check()

                parsed_value = normpath(value) if dtype == 'path' else dtype(
                    value)
                config_args.update({key: parsed_value})

        return config_args

    def check_constraints(self):
        raise NotImplementedError


class STDConfig(Config):

    def check_constraints(self):
        self.build_dir_hierarchy()
        self.config_args['std-spectra'] = self.check_spectra()
        self.config_args['std-regions'] = self.check_regions()

    def build_dir_hierarchy(self):
        # Todo: exist_ok=False for all occurrences
        os.makedirs(self.config_args['output_dir'], exist_ok=True)

        outdir = self.config_args['output_dir']
        for dir_name in ['STD', 'DOCKING']:
            self.config_args[dir_name] = join(outdir, dir_name)
            os.makedirs(self.config_args[dir_name], exist_ok=True)

        # Write the configuration file for reproducibility
        stdock_config = join(self.config_args['output_dir'], 'stdock.config')
        with open(stdock_config, 'wt') as ini:
            self.config_obj.write(ini)

    def check_spectra(self):
        spectra = list(self.config_obj['std-spectra'].keys())
        spectra_dict = {}
        for i, pair in enumerate(spectra):
            on, off, d20 = map(str.strip, pair.split(','))
            d20 = float(d20)
            cmn.check_numeric_in_range('d20', d20, float, 0, np.inf)
            on_present = os.path.exists(on)
            off_present = os.path.exists(off)

            if not (on_present and off_present):
                raise ValueError(f'Either {on} or {off} is not a valid path')
            spectra_dict.update({i: {'on': on, 'off': off, 'd20': d20}})
        return spectra_dict

    def check_regions(self):
        regions = list(self.config_obj['std-regions'].keys())
        regions_dict = {}
        for i, entry in enumerate(regions):
            label, upper_raw, lower_raw = entry.split(',')
            regions_dict.update({i: {'label': label, 'upper': float(upper_raw),
                                     'lower': float(lower_raw)}})
        return regions_dict


# =============================================================================
# Debugging area
# =============================================================================
# config_path = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/config.cfg'
# params = allowed_parameters
# templates = allowed_templates
# self = STDConfig(config_path, params, templates)
