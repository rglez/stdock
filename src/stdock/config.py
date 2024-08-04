# Created by roy.gonzalez-aleman at 13/11/2023
import configparser
import os
import textwrap
from collections import Counter
from os.path import abspath, dirname, isabs, join, normpath

import stdock.commons as cmn

#: Total number of cores
num_cores = os.cpu_count()

#: Allowed section templates in the config file
allowed_templates = {
    'map-from-spectra':
        {'generals', 'std-regions', 'std-spectra'},

    'map-from-values-then-dock-small':
        {'generals', 'std-epitope', 'docking-small'},

    'map-from-spectra-then-dock-small':
        {'generals', 'std-regions', 'std-spectra', 'docking-small'},

    'map-from-spectra-then-dock-pept':
        {'generals', 'std-regions', 'std-spectra', 'docking-pept'},

    # 'dock-without-map':
    #     {}
}

#: Allowed keys in the config file (dtypes & expected values)
allowed_parameters = {

    # General parameters
    'generals': {
        'output_dir': {'dtype': 'path', 'check_exist': False},
        'protein_conc': {'dtype': float, 'min': 0.0, 'max': cmn.inf_float}},

    # STD-NMR integral regions parameters
    'std-regions': None,

    # (On-res, off-res, spectra) tuples
    'std-spectra': None,

    # STD-NMR epitope mapping parameters
    'std-epitope': None,

    # VINADocking-related parameters
    'docking-small': {
        'receptor_pdb': {'dtype': 'path', 'check_exist': True},
        'ligand_pdb': {'dtype': 'path', 'check_exist': True},
        'num_poses': {'dtype': int, 'min': 1, 'max': cmn.inf_int},
        'rmsd_tolerance': {'dtype': float, 'min': 0.01, 'max': cmn.inf_float},
        'exhaustiveness': {'dtype': int, 'min': 1, 'max': cmn.inf_int}},

    'docking-pept': {
        'receptor_pdb': {'dtype': 'path', 'check_exist': True},
        'ligand_pdb': {'dtype': 'path', 'check_exist': True},
        'scoring_function': {'dtype': str,
                             'values': {'cpydock', 'dfire', 'fastdfire',
                                        'dfire2', 'dna', 'ddna',
                                        'mj3h', 'pisa', 'sd', 'sipper', 'tobi',
                                        'vdw'}},
        'num_steps': {'dtype': int, 'min': 1, 'max': cmn.inf_int}
    }
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
        self.parse_and_check_constraints()
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
        [current_template.remove(x) for x in self.keyless_sections if
         x in current_template]

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
        [parsed_sections.remove(x) for x in self.keyless_sections if
         x in parsed_sections]

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

    def parse_and_check_constraints(self):
        """Check for constraints in the STDock config file
        """
        raise NotImplementedError


class STDConfig(Config):
    """
    Specific parser for STDock's config files. It inherits from a more general
    config parser and then perform STDock-related checkings.
    """

    def parse_and_check_constraints(self):
        config_sections = self.config_obj.sections()
        # 1. Build dir hierarchy
        self.build_dir_hierarchy()

        # 2. Check [std-spectra] section
        if 'std-spectra' in config_sections:
            self.config_args['std-spectra'] = self.parse_spectra()

        # 3. Check [std-regions] section
        if 'std-regions' in config_sections:
            self.config_args['std-regions'] = self.parse_regions()

        # 4. Check [std-epitopes] section
        if 'std-epitope' in config_sections:
            self.config_args['std-epitopes'] = self.parse_epitopes()

        # 5. Get protein concentration in general section if specified
        # prot_conc = self.config_obj['generals'].get('protein_conc', None)
        # self.config_args['prot_conc'] = float(prot_conc) if prot_conc else None

        # todo: Check [docking] section

    def build_dir_hierarchy(self):
        """
        Build STDock directory hierarchy
        """
        # If output_dir exists, raise
        outdir = self.config_args['output_dir']
        try:
            os.makedirs(outdir)
        except FileExistsError:
            raise FileExistsError(
                f'The output directory {outdir} already exists. Please, '
                f'choose another one, or delete the existing one.')
        for dir_name in ['STD', 'DOCKING']:
            self.config_args[dir_name] = join(outdir, dir_name)
            os.makedirs(self.config_args[dir_name])

        # Write the configuration file for reproducibility
        stdock_config = join(self.config_args['output_dir'], 'stdock-job.cfg')
        with open(stdock_config, 'wt') as ini:
            self.config_obj.write(ini)

    def parse_spectra(self):
        """
        Parse the [std-spectra] section

        Returns:
            spectra_dict: a nested dict with ligand concentrations and
                          saturation times respectively
        """
        # Extract keys and values
        config = self.config_obj
        keys = list(config['std-spectra'].keys())
        values_raw = [x.strip() for x in list(config['std-spectra'].values())]

        # Perform path validity checks for spectra files
        values = [cmn.check_path(x) for x in values_raw]

        # Get types, times & concentrations in the [std-spectra] section
        types, times, concs = [], [], []
        for key in keys:
            splitted = key.split('_')
            if len(splitted) != 3:
                raise ValueError(
                    'Keys of the [std-spectra] section of the config file must'
                    ' be labeled as [on, off or diff]_[ligand concentration]_[saturation time].'
                    ' For example: on_0.1_0.5, off_0.1_0.5, diff_0.1_0.5')

            types.append(splitted[0])
            times.append(float(splitted[1]))
            concs.append(float(splitted[2]))

        # Check types of spectra are: on, off or diff
        types_counts = Counter(types)
        valid_types = ['on', 'off', 'diff']
        types_are_valid = all([x in valid_types for x in types_counts.keys()])
        if not types_are_valid:
            raise ValueError(
                f'Invalid spectrum types in the [std-spectra] of the '
                f'config file. Valid types are: {valid_types}')

        # Check only two types are present
        if len(types_counts) != 2:
            raise ValueError(
                'There must be only two types of spectrum specified in the '
                '[std-spectra] section of the config file: (on and off) OR'
                ' (diff and off)')

        # Check types number consistency
        if len(set(types_counts.values())) != 1:
            raise ValueError(
                'There must be the same number of (on and off) OR (diff and off)'
                ' spectra specified in the [std-spectra] section of the config'
                ' file')

        # Check types consistency
        conc_times = Counter(list(zip(concs, times))).values()
        if set(conc_times) != {2}:
            raise ValueError(
                'There must be two spectrum per each'
                ' (concentration/saturation time) level')

        # Check times consistency
        if len(set(times)) < 2:
            raise ValueError(
                'There must be at least two different saturation times for each '
                'ligand concentration specified')

        # Check concs consistency
        if len(set(concs)) < 1:
            raise ValueError(
                'There must be at least one ligand concentration specified')

        if len(set(Counter(concs).values())) != 1:
            raise ValueError(
                'There must be the same number of spectra for each ligand '
                'concentration specified in the [std-spectra] section of the '
                'config file')

        # Get all the spectra organized as a dict
        spectra_dict = cmn.recursive_defaultdict()
        for i, key in enumerate(keys):
            splitted = key.split('_')
            type_ = splitted[0]
            lig_conc, sat_time = map(float, splitted[2:0:-1])
            spectra_dict[lig_conc][sat_time][type_] = values[i]
        return spectra_dict

    def check_spectra(self):
        # todo: check all relative to the spectra consistency
        pass

    def parse_regions(self):
        """
        Parse the [std-regions] section

        Returns:
            regions_dict: a dict with regions limits to integrate
        """
        error_regions = textwrap.fill(
            '''
            Integration regions in the config file must be specified as {x}_{y}_{z}
            one-liners where x is an atomic label and y, z are the upper and lower
            limits of the integration region respectively''', width=80)

        # Extract keys and values
        config = self.config_obj
        keys = list(config['std-regions'].keys())
        values_raw = list(config['std-regions'].values())

        # Get all the integral regions organized as a dict
        regions_dict = {}
        for i, key in enumerate(keys):
            upper_raw, lower_raw = values_raw[i].split(',')
            if (upper := float(upper_raw)) > (lower := float(lower_raw)):
                regions_dict.update(
                    {i: {'label': key, 'upper': upper, 'lower': lower}})
            else:
                raise ValueError(error_regions)
        return regions_dict

    def parse_epitopes(self):
        config = self.config_obj
        keys = list(config['std-epitope'].keys())
        values_raw = [float(x) for x in config['std-epitope'].values()]
        epitope = dict(zip(keys, values_raw))
        return epitope

# %%===========================================================================
# Debugging area
# =============================================================================
# Debugging 'map-from-spectra'
# config_path = '/home/gonzalezroy/RoyHub/stdock/tests/paper/imaging/config_kd.cfg'
# params = allowed_parameters
# templates = allowed_templates
# self = STDConfig(config_path, params, templates)
# self.config_args

# Debugging 'map-from-spectra-then-dock'
# config_path = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/troll8.cfg'
# params = allowed_parameters
# templates = allowed_templates
# self = STDConfig(config_path, params, templates)

# Debugging 'map-from-values-then-dock'
# config_path = '/home/roy.gonzalez-aleman/RoyHub/stdock/tests/example/00-CASE_STUDY/HuR/M9/M9.cfg'
# params = allowed_parameters
# templates = allowed_templates
# self = STDConfig(config_path, params, templates)
