# Created by roy.gonzalez-aleman at 13/11/2023
import fnmatch
import os
import pickle
import sys
from collections import defaultdict

import numpy as np


def check_path(path, check_exist=True):
    path_exists = os.path.exists(path)
    if check_exist and path_exists:
        return path
    elif (not check_exist) and (not path_exists):
        return path
    elif (not check_exist) and path_exists:
        return path
    elif check_exist and (not path_exists):
        raise ValueError(f'\nNo such file or directory: {path}')
    else:
        pass
        raise ValueError(
            f'\nPath already exists and will not be overwritten: {path}')


def recursive_finder(pattern, root=os.curdir):
    for path, dirs, files in os.walk(os.path.abspath(root), followlinks=True):
        if dirs:
            for dir_ in dirs:
                recursive_finder(dir_)
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(path, filename)


def recursive_defaultdict():
    return defaultdict(recursive_defaultdict)


def check_file_extension(path, extension):
    ext = path.split(os.sep)[-1].split('.')[-1]
    if ext != extension:
        raise TypeError(f'\n{path} file extension is not {extension}')


def check_numeric_in_range(arg_name, value, dtype, minim, maxim):
    if not isinstance(value, dtype):
        raise TypeError(f'Param "{arg_name}" must be of type {dtype}')

    if not minim <= value <= maxim:
        raise ValueError(f'Param "{value}" out of [{minim}, {maxim}]')

    return dtype(value)


def combinator(top_iterable):
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
    ''' Serialize data using **pickle**.

    Args:
        data (object)  : any serializable object.
        file_name (str): name of the **pickle** file to be created.
    Returns:
        (str): file_name
    '''
    with open(file_name, 'wb') as file:
        pickle.dump(data, file)
    return file_name


def unpickle_from_file(file_name):
    ''' Unserialize a **pickle** file.

    Args:
        file_name (str): file to unserialize.
    Returns:
        (object): an unserialized object.
    '''
    with open(file_name, 'rb') as file:
        data = pickle.load(file)
    return data


inf_int = sys.maxsize
inf_float = float(inf_int)
