# Created by roy.gonzalez-aleman at 13/11/2023
import fnmatch
import os
import pickle
import shutil
import subprocess as sp
import sys
import tempfile
from collections import defaultdict
from os.path import join, split

import matplotlib as mpl

mpl.use('TKAgg')

import mdtraj as md
import numpy as np
import prody as prd
from numba import jit, prange
from openbabel import pybel
from scipy.spatial import cKDTree as ckd

ob_log_handler = pybel.ob.OBMessageHandler()
ob_log_handler.SetOutputLevel(0)

inf_int = sys.maxsize
inf_float = float(inf_int)

syn = {8: 'low',
       64: 'medium',
       256: 'high',
       'speed1': 'high',
       'speed2': 'medium',
       'speed4': 'low',
       125000: 'low',
       250000: 'medium',
       500000: 'high'}


def check_epitope_label(label):
    """
    Check if a string is a valid label

    Args:
        label: a string to check

    Returns:
        label: the checked label
    """
    try:
        atomic_label, user_label = label.split('-')
        atom_name, atom_number = atomic_label.split('_')
    except ValueError:
        atomic_label, user_label = None, None
        raise ValueError(
            'Labels of sections [std-regions] and [std-epitope] must have the'
            f' format: [atom_name]_[atom_number]-[user_label].'
            f' You passed {label}')
    return label


def check_path(path, check_exist=True):
    """
    Check if a path exists and return it if it does

    Args:
        path: path to check
        check_exist: check if the path exists or not

    Returns:
        path: path to the file or directory
    """
    path_exists = os.path.exists(path)
    if check_exist and path_exists:
        return path
    elif (not check_exist) and (not path_exists):
        return path
    elif (not check_exist) and path_exists:
        return path  # todo: check this behaviour
    elif check_exist and (not path_exists):
        raise ValueError(f'\nNo such file or directory: {path}')
    else:
        pass
        raise ValueError(
            f'\nPath already exists and will not be overwritten: {path}')


def pdb2pdbqt(pdb_path, ligand_case=False,
              babel_path=shutil.which('obabel')):
    """
    Convert a pdb file to pdbqt format using OpenBabel

    Args:
        pdb_path: path to the pdb file
        ligand_case: is the pdb a ligand or not?
        babel_path: path to the obabel executable

    Returns:
        pdbqt_path: path to the pdbqt file
    """
    dirname, basename_raw = split(pdb_path)
    basename = join(dirname, basename_raw.split('.')[0])

    if not ligand_case:
        cmd = (f'{babel_path} -ipdb {pdb_path} -opdbqt -xp -xh -xn -xr -O'
               f' {basename}.pdbqt --partialcharge gasteiger')
    else:
        cmd = (f'{babel_path} -ipdb {pdb_path} -opdbqt -xp -xh -xn -O'
               f' {basename}.pdbqt --partialcharge gasteiger')
    cmd_list = cmd.split()
    output, errors = shell_run(cmd_list)
    return check_path(f'{basename}.pdbqt')


def shell_run(cmd_list):
    """
    Run a command in a subprocess

    Args:
        cmd_list: command to run split as list

    Returns:
        outputs: stdout
        errors: stderr
    """
    cmd_run = sp.Popen(cmd_list, text=False, stdout=sp.PIPE, stderr=sp.PIPE)
    output, errors = cmd_run.communicate()
    return output, errors


def makedir_after_overwriting(path):
    """
    Make a directory after remooving it if already created

    Args:
        path: path to the dir to be (re)-created

    Returns:
        path: path to the dir
    """
    shutil.rmtree(path, ignore_errors=True)
    os.makedirs(path)
    return path


def recursive_finder(pattern, root=os.curdir):
    """
    Find all files in a directory tree that match a pattern

    Args:
        pattern: pattern to match
        root: root directory to start the search

    Returns:
        iterator: all files that match the pattern
    """
    for path, dirs, files in os.walk(os.path.abspath(root), followlinks=True):
        if dirs:
            for dir_ in dirs:
                recursive_finder(dir_)
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(path, filename)


def recursive_defaultdict():
    """
    Create a recursive defaultdict

    Returns:
        defaultdict: recursive defaultdict
    """
    return defaultdict(recursive_defaultdict)


# def check_file_extension(path, extension):
#     ext = path.split(os.sep)[-1].split('.')[-1]
#     if ext != extension:
#         raise TypeError(f'\n{path} file extension is not {extension}')


def check_numeric_in_range(arg_name, value, dtype, minim, maxim):
    """
    Check if a value is of a certain type and within a range

    Args:
        arg_name: name of the argument to check
        value: value to check
        dtype: type of the value
        minim: minimum value
        maxim: maximum value

    Returns:
        value: value of the correct type and within the range
    """
    if not isinstance(value, dtype):
        raise TypeError(f'Param "{arg_name}" must be of type {dtype}')

    if not minim <= value <= maxim:
        raise ValueError(f'Param "{value}" out of [{minim}, {maxim}]')

    return dtype(value)


def combinator(top_iterable):
    """
    Create all possible combinations of a list of iterables

    Args:
        top_iterable: list of iterables

    Returns:
        ANG: list of all possible combinations
    """
    #  Constructors parameters
    LEN = [len(x) for x in top_iterable]
    N = np.prod(LEN)
    DIV = [np.prod(LEN[0:x + 1]) for x in range(len(top_iterable))]
    REP = [int(N / D) for D in DIV]
    repetitions = [int(N / (REP[x] * LEN[x])) for x in
                   range(len(top_iterable))]
    # Combinator
    columns = []
    for index, iterable in enumerate(top_iterable):
        col = []
        for idx, element in enumerate(iterable):
            r = REP[index]
            while r != 0:
                col.append(element)
                r -= 1
        columns.append(col)

    #  Final product
    COMB = [iterable * repetitions[index]
            for index, iterable in enumerate(columns)]

    #  List of angles creation
    ANG = []
    for index in range(len(COMB[0])):
        conformer = []
        for idx, iterable in enumerate(COMB):
            conformer.append(COMB[idx][index])
        ANG.append(conformer)
    return ANG


def pickle_to_file(data, file_name):
    """
    Save a data structure to a file using pickle

    Args:
        data: data to save
        file_name: path to the file to save the data

    Returns:
        file_name: path to the file
    """
    with open(file_name, 'wb') as file:
        pickle.dump(data, file)
    return file_name


def unpickle_from_file(file_name):
    """
    Load a data structure from a file using pickle

    Args:
        file_name: path to the file to load the data

    Returns:
        data: data loaded from the file
    """
    with open(file_name, 'rb') as file:
        data = pickle.load(file)
    return data


def generic_matplotlib():
    """
    Set generic values for matplotlib-generated plots
    """
    mpl.rc('figure', figsize=[12, 8], dpi=600)
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
    """
    Reset matplotlib values to defaults
    """
    mpl.rcParams.update(mpl.rcParamsDefault)


def get_final_names(mapping_labels):
    """
    Get the final names of the columns in the dataframe

    Args:
        mapping_labels: list of column names

    Returns:
        final_names: list of final column names
    """
    exh_dict = {'high': 'H', 'medium': 'M', 'low': 'L'}
    sf_dict = {'dkoes_fast': 'dkoes1',
               'dkoes_scoring': 'dkoes2',
               'ad4_scoring': 'ad4',
               'autodock4': 'ad4'}
    prog_dict = {'autodock4': 'ad4'}
    final_names = []
    for i, name in enumerate(mapping_labels):
        splitted = name.split('_')
        N = len(splitted)
        if N == 3:
            prog, sf, exh = splitted
        elif N == 4:
            prog, *sf, exh = splitted
            sf = '_'.join(sf)
        else:
            raise ValueError('Error splitting column name')
        prog = prog_dict.get(prog, prog)
        sf = sf_dict.get(sf, sf)
        exh = exh_dict[exh]
        final_names.append((prog.upper(), sf, exh))
    return final_names


# def chop_cmap_frac(cmap, frac):
#     """Chops off the beginning `frac` fraction of a colormap."""
#     cmap_as_array = cmap(np.arange(256))
#     cmap_as_array = cmap_as_array[int(frac * len(cmap_as_array)):]
#     list_name = cmap.name + f"_frac{frac}"
#     return LinearSegmentedColormap.from_list(list_name, cmap_as_array)


def get_filtered_indices(rec_kdt, lig_parsed):
    """
    Filter ligands by distance to the receptor

    Args:
        rec_kdt: cKDTree object of the receptor
        lig_parsed: parsed ligands

    Returns:
        filtered_indices: indices of the ligands that are close to the receptor
        filtered_ligs: ligands that are close to the receptor
    """
    filtered_indices = []
    filtered_ligs = []
    for i, lig_ag in enumerate(lig_parsed):
        try:
            lig_kdt = ckd(lig_ag.getCoords())
            contacts = lig_kdt.query_ball_tree(rec_kdt, r=5)
            if any(contacts):
                filtered_indices.append(i)
                filtered_ligs.append(lig_ag)
        except AttributeError:
            print(f'WARNING: Frame {i} not filtered')
    return filtered_indices, filtered_ligs


def write_string(string, path):
    """
    Write a string to a file

    Args:
        string: string to write
        path: path to the file to write the string

    Returns:
        path: path to the file
    """
    with open(path, 'wt') as out:
        out.write(string)
    return path


def get_longest_component(parsed_mol):
    """
    Get the size of the longest component of a molecule

    Args:
        parsed_mol: parsed molecule

    Returns:
        size: size of the longest component
    """
    coords = parsed_mol.getCoords()
    mini = coords.min(axis=0)
    maxi = coords.max(axis=0)
    size = abs((maxi - mini).max())
    return size


def reformat_single_x2y(x_mol_path, y_mol_path):
    """
    Reformat a molecule from one format

    Args:
        x_mol_path: path to the molecule in format x
        y_mol_path: path to the molecule in format y
    """
    x_root, x_base = split(x_mol_path)
    y_root, y_base = split(y_mol_path)

    x_format = x_base.split('.')[-1]
    y_format = y_base.split('.')[-1]

    parsed = next(pybel.readfile(x_format, x_mol_path))

    parsed.write(y_format, y_mol_path)


class Molecule:
    """
    Class to parse a molecule from a path
    """

    def __init__(self, molecule_path):
        self.mol_path = molecule_path
        check_path(self.mol_path)

    def parse(self):
        """
        Parse a molecule from a path

        Returns:
            parsed: parsed molecule
        """
        root, base = split(self.mol_path)
        extension = base.split('.')[-1]
        rec_obj = pybel.readfile(extension, self.mol_path)

        with tempfile.TemporaryDirectory() as d:
            parsed = []
            for i, obj in enumerate(rec_obj):
                temp_mol = join(d, f'temp_molecule_{i}.pdb')
                obj.write("pdb", temp_mol)
                parsed.append(prd.parsePDB(temp_mol))
        return parsed

    def get_ensemble(self):
        """
        Get a Prody ensemble of the parsed molecule

        Returns:
            ensemble: Prody ensemble of the parsed molecule
        """
        parsed = self.parse()
        ensemble = prd.Ensemble()
        for i, ag in enumerate(parsed):
            try:
                ensemble.addCoordset(ag.getCoords())
            except AttributeError:
                print(f'WARNING: Frame {i} not parsed.')
                pass
        ensemble.setAtoms(parsed[0])
        return ensemble


class Program:
    """
    Prototype of a Program class
    """

    def __init__(self, exe_path, rec_path, lig_path, n_poses, rmsd_tol,
                 scoring_functions, exhaustiveness_list, odir):
        # Program executable path
        self.exe = check_path(exe_path)

        # Receptor & Ligand parsing
        self.rec_path = check_path(rec_path)
        self.rec_parsed = Molecule(self.rec_path).parse()[0]
        self.lig_path = check_path(lig_path)
        self.lig_parsed = Molecule(self.lig_path).parse()[0]

        # Gather all scoring functions
        self.scoring_functions = scoring_functions

        # Exhaustiveness
        self.exhaustiveness_list = exhaustiveness_list

        # Check numeric arguments
        self.num_poses = check_numeric_in_range('n_poses', n_poses, int, 1,
                                                inf_int)
        self.rmsd_tol = check_numeric_in_range('rmsd_tol', rmsd_tol, float,
                                               0.1, inf_float)
        # Create output dir
        self.out_dir = odir
        os.makedirs(self.out_dir, exist_ok=True)

        # Standardize programs' running
        self.commands, self.config_paths = self.get_commands()
        self.bench = ['/usr/bin/time', '-v']

    def get_rec_axis(self, output_min_max=False):
        """
        Get the axis of the receptor

        Args:
            output_min_max: output the min and max values of the axis?

        Returns:
            axis: axis of the receptor
        """
        rec_coords = self.rec_parsed.getCoords()
        min_coords = rec_coords.min(axis=0)
        max_coords = rec_coords.max(axis=0)
        if output_min_max:
            return max_coords - min_coords, min_coords, max_coords
        return max_coords - min_coords

    def get_rec_center(self):
        """
        Get the center of the receptor

        Returns:
            center: center of the receptor

        """
        rec_parsed = Molecule(self.rec_path).parse()[0]
        center = np.round(prd.calcCenter(rec_parsed), 1)
        return center

    def get_commands(self):
        """
        Get the commands to run the program
        """
        raise NotImplementedError

    def run_commands(self):
        """
        Run the commands to dock the ligand
        """
        raise NotImplementedError

    def yield_filter_sort(self):
        """
        Yield the filtered and sorted results
        """
        raise NotImplementedError


@jit(nopython=True, parallel=True)
def rmsd_ref_vs_all(ref, array, div_by_n):
    """
    Calculate the RMSD between a reference and an array of coordinates

    Args:
        ref: reference coordinates
        array: array of coordinates
        div_by_n: 1/N value

    Returns:
        rmsds: array of RMSD values
    """
    rmsds = np.zeros(len(array), dtype=float)
    for i in prange(len(array)):
        tar = array[i]
        rmsds[i] = np.sqrt(((ref - tar) ** 2).sum() * div_by_n)
    return rmsds


def leader_clustering_matrix(sorted_ensemble, cutoff):
    """
    Cluster a matrix using the leader algorithm

    Args:
        sorted_ensemble: ensemble of conformations
        cutoff: RMSD cutoff

    Returns:
        clusters: list of clusters
    """
    clusters = []
    N = len(sorted_ensemble)
    clustered = np.zeros(N, dtype=bool)
    array = sorted_ensemble.getCoordsets()
    divByN = 1.0 / array[0].shape[0]

    while not clustered.all():
        next_frame = clustered.argmin()
        rmsd = rmsd_ref_vs_all(array[next_frame], array, divByN)
        new_cluster = rmsd <= cutoff
        true_clusters = np.bitwise_and(new_cluster, ~clustered)
        clustered[new_cluster] = True
        clusters.append(true_clusters.nonzero()[0])
    return clusters


def leader_clustering_traj(sorted_dcd_path, sorted_pdb_path, cutoff):
    """
    Cluster a trajectory using the leader algorithm

    Args:
        sorted_dcd_path: path to the sorted dcd file
        sorted_pdb_path: path to the sorted pdb file
        cutoff: rmsd cutoff

    Returns:
        clusters: list of clusters
    """
    if sorted_dcd_path.split('.')[-1] == 'dcd':
        trajectory = md.load(sorted_dcd_path, top=sorted_pdb_path)
    elif sorted_dcd_path.split('.')[-1] == 'pdb':
        trajectory = md.load(sorted_dcd_path)
    else:
        raise ValueError(f'Only pdb or dcd formats are available for traj')
    trajectory.center_coordinates()
    clusters = []
    N = trajectory.n_frames
    clustered = np.zeros(N, dtype=bool)

    while not clustered.all():
        next_frame = clustered.argmin()
        rmsd = md.rmsd(trajectory, trajectory, next_frame, precentered=True)
        new_cluster = rmsd <= cutoff
        true_clusters = np.bitwise_and(new_cluster, ~clustered)
        clustered[new_cluster] = True
        clusters.append(true_clusters.nonzero()[0])
    return clusters
