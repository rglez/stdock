# Created by roy.gonzalez-aleman at 03/12/2023
import argparse
import os
import pickle
from collections import deque

import numpy as np
import pandas as pd
import prody as prd
from bitarray import bitarray as ba
from bitarray import util as bu


def binarize_matrix(arguments, float_matrix):
    N = float_matrix.shape[0]
    cutoff = np.full(N, arguments.cutoff, dtype=np.float32)
    matrix = dict()
    for i in range(N):
        rmsd_ = float_matrix[i]
        vector_np = np.less_equal(rmsd_, cutoff)
        bit_array = ba()
        bit_array.pack(vector_np.tobytes())
        bit_array.fill()
        matrix.update({i: bit_array})
    return matrix


def calc_matrix_degrees(unclustered_bit, matrix):
    """
    Calculate number of neighbors (degree) of unclustered nodes in matrix.

    Parameters
    ----------
    unclustered_bit : bitarray.bitarray
        bitarray with indices of unclustered nodes turned on.
    matrix : collections.OrderedDict
        dict of bitarrays.

    Returns
    -------
    degrees : numpy.ndarray
        array containing each node degree. Clustered nodes have degree = 0.

    """
    one = ba('1')
    degrees = np.zeros(len(unclustered_bit), dtype=np.int32)
    for node in unclustered_bit.itersearch(one):
        try:
            degrees[node] = matrix[node].count()
        except KeyError:
            pass
    return degrees


def colour_matrix(degrees, matrix):
    """
    Greedy coloring of bit-encoded RMSD matrix.

    Parameters
    ----------
    degrees : numpy.ndarray
        array containing each node degree. Clustered nodes have degree = 0.
    matrix : collections.OrderedDict
        dict of bitarrays.

    Returns
    -------
    colors : numpy.ndarray
        array of colors assigned to each node of the matrix.
    """
    # Constants ---------------------------------------------------------------
    N = degrees.size
    m = len(matrix)
    one = ba('1')
    xcolor = 0
    # Initialize containers ---------------------------------------------------
    ordered_by_degrees = iter((-degrees[:m]).argsort())
    colors = np.zeros(N, dtype=np.int32)
    colored = ba(N)
    colored.setall(0)
    seen = set()
    while True:
        # Retrieve the max-degree node ----------------------------------------
        max_node = next(ordered_by_degrees)
        if max_node in seen:
            continue
        seen.add(max_node)
        xcolor += 1
        not_neighbors = ~ matrix[max_node]
        not_colored = ~colored
        candidates = not_neighbors & not_colored
        # Nodes passing conditions (not-neighb, not-colored, not-neighb) ------
        passed = [max_node]
        for candidate in candidates.itersearch(one):
            passed.append(candidate)
            try:
                candidates &= ~matrix[candidate]
            except KeyError:
                continue
            if not candidates.any():
                break
        seen.update(passed)
        # Deliver a color class to passed nodes -------------------------------
        colors[passed] = xcolor
        colored = ba()
        colored.pack(colors.astype(np.bool).tobytes())
        if colored.count(0) == 0:
            break
    return colors


def bitarray_to_np(bitarr):
    """
    Convert from bitarray.bitarray to numpy.ndarray efficiently.

    Parameters
    ----------
    bitarr : bitarray.bitarray
        a bitarray.

    Returns
    -------
    numpy.ndarray
        boolean bitarray equivalent to the binary bitarray input object.
    """
    return np.unpackbits(bitarr).astype(bool)


def do_bit_cascade(big_node, degrees, colors, matrix, max_):
    """
    Perform succesive AND operations between an initial bitarray and subsequent
    bitarray candidates to search for a clique.

    Parameters
    ----------
    big_node : int
        node whose bitarray will start the operations.
    degrees : numpy.ndarray
        array containing each node degree. Clustered nodes have degree = 0.
    colors : numpy.ndarray
        array of colors assigned to each node of the matrix.
    clustered_bit : bitarray.bitarray
        bitarray with indices of clustered nodes turned on.
    matrix : collections.OrderedDict
        dict of bitarrays.
    max_ : int
        Stop iterative AND operations after the initial bitarray has max_
        bits turned on.

    Returns
    -------
    init_cascade : bitarray.bitarray
        initial bitarray before any AND operation.
    ar : numpy.ndarray
        array of nodes forming a clique.
    """
    init_cascade = matrix[big_node]
    # .... recovering neighbors and their information .........................
    neighb = bitarray_to_np(init_cascade).nonzero()[0]
    neighb_colors = colors[neighb]
    if len(set(neighb_colors.tolist())) <= max_:
        return None
    neighb_degrees = degrees[neighb]
    g = np.bincount(neighb_colors)
    neighb_g = g[neighb_colors]
    # .... ordering neighbors by g ---> colors ---> degrees ...................
    idx = np.lexsort([-neighb_degrees, neighb_colors, neighb_g])
    candidates_info = zip(neighb[idx], neighb_colors[idx])

    # .... BitCascade considering divergence ..................................
    counter = 0
    seen = set()
    for candidate, color in candidates_info:
        if (color in seen) or (not init_cascade[candidate]):
            continue
        seen.add(color)
        init_cascade = matrix[candidate] & init_cascade
        counter += 1
        COUNT = init_cascade.count()
        if (COUNT <= max_):
            return None
        if counter >= COUNT:
            break
    ar = np.nonzero(np.unpackbits(init_cascade).astype(bool))[0]
    return init_cascade, ar


def set_to_bitarray(set_, N):
    """
    Convert from python set to bitarray.bitarray.

    Parameters
    ----------
    set_ : set
        a python set.
    N : int
        lenght of the desired bitarray. It must be greater than the maximum
        value of indices present in set.

    Returns
    -------
    bitarr : bitarray.bitarray
        bitarray of lenght N with indices present in set turned on.
    """
    zero_arr = np.zeros(N, dtype=bool)
    zero_arr[list(set_)] = 1
    bitarr = ba()
    bitarr.pack(zero_arr.tobytes())
    return bitarr


def pickle_to_file(data, file_name):
    """
    Serialize data using the pickle library.

    Parameters
    ----------
    data : serializable object
        variable name of the object to serialize.
    file_name : str
        name of the pickle file to be created.

    Returns
    -------
    file_name : str
        name of the serialized file.
    """
    with open(file_name, 'wb') as file:
        pickle.dump(data, file)
    return file_name


def top_has_coords(topology):
    """
    Check if topology has cartesian coordinates information.

    Parameters
    ----------
    topology : str
        Path to the topology file.

    Returns
    -------
    int
        Number of cartesian frames if topology contains cartesians.
        False otherwise.
    """
    try:
        tt = md.load(topology)
    except OSError:
        return False
    return tt.xyz.shape[0]


def to_VMD(outdir, topology, first, N1, last, stride, final_array):
    """
    Create a .log file for visualization of clusters in VMD through a
    third-party plugin.

    Parameters
    ----------
    outdir : str
        Path where to create the VMD visualization .log.
    topology : str
        Path to the topology file.
    first : int
        First frame to consider (0-based indexing).
    N1 : int
        default value when last == None.
    last : TYPE
        Last frame to consider (0-based indexing).
    stride : TYPE
        Stride (step).
    final_array : numpy.ndarray
        Final labeling of the selected clusters ordered by size (descending).

    Returns
    -------
    logname : str
        Log file to be used with VMD.
    """
    basename = os.path.basename(topology).split('.')[0]
    logname = os.path.join(outdir, '{}.log'.format(basename))
    vmd_offset = top_has_coords(topology)
    start = first
    if not last:
        stop = N1
    else:
        stop = last
    slice_frames = np.arange(start, stop, stride, dtype=np.int32)
    nmr_offset = 1
    with open(logname, 'wt') as clq:
        for num in np.unique(final_array):
            if num != 0:
                clq.write('{}:\n'.format(num))
                cframes = np.where(final_array == num)[0]
                if vmd_offset:
                    real_frames = slice_frames[
                                      cframes] + nmr_offset + vmd_offset
                else:
                    real_frames = slice_frames[cframes] + nmr_offset
                str_frames = [str(x) for x in real_frames]
                members = ' '.join(str_frames)
                clq.write('Members: ' + members + '\n\n')
        if 0 in np.unique(final_array):
            clq.write('{}:\n'.format(0))
            cframes = np.where(final_array == 0)[0]
            if vmd_offset:
                real_frames = slice_frames[cframes] + nmr_offset + vmd_offset
            else:
                real_frames = slice_frames[cframes] + nmr_offset
            str_frames = [str(x) for x in real_frames]
            members = ' '.join(str_frames)
            clq.write('Members: ' + members + '\n\n')
    return logname


def get_frames_stats(N1, first, last, stride, clusters, outdir):
    """
    Get "frames_statistics.txt" containing frameID, clusterID.

    Parameters
    ----------
    N1 : int
        default value when last == None.
    first : int
        First frame to consider (0-based indexing).

    last : TYPE
        Last frame to consider (0-based indexing).
    stride : TYPE
        Stride (step).
    clusters : numpy.ndarray
        array of clusters ID.
    outdir : str
        Path where to create the VMD visualization .log.

    Returns
    -------
    frames_df : pandas.DataFrame
        dataframe with frames_statistics info.
    """
    start = first
    if not last:
        stop = N1
    else:
        stop = last
    slice_frames = np.arange(start, stop, stride, dtype=np.int32)
    frames_df = pd.DataFrame(columns=['frame', 'cluster_id'])
    frames_df['frame'] = range(N1)
    frames_df['cluster_id'].loc[slice_frames] = clusters
    with open(os.path.join(outdir, 'frames_statistics.txt'), 'wt') as on:
        frames_df.to_string(buf=on, index=False)
    return frames_df


def get_cluster_stats(clusters, outdir):
    """
    Get "cluster_statistics.txt" containing clusterID, cluster_size, and
    cluster percentage from trajectory.

    Parameters
    ----------
    clusters : numpy.ndarray
        array of clusters ID.
    outdir : str
        Path where to create the VMD visualization .log.

    Returns
    -------
    clusters_df : pandas.DataFrame
        dataframe with cluster_statistics info.
    """
    clusters_df = pd.DataFrame(columns=['cluster_id', 'size', 'percent'])
    clusters_df['cluster_id'] = list(range(0, clusters.max() + 1))
    sizes = []
    for x in clusters_df.cluster_id:
        sizes.append(len(np.where(clusters == x)[0]))
    clusters_df['size'] = sizes

    sum_ = clusters_df['size'].sum()
    percents = [round(x / sum_ * 100, 4) for x in clusters_df['size']]
    clusters_df['percent'] = percents

    with open(os.path.join(outdir, 'cluster_statistics.txt'), 'wt') as on:
        clusters_df.to_string(buf=on, index=False)
    return clusters_df


def bitqt_matrix(args, float_matrix):
    # =========================================================================
    # 1. Creating binary matrix (adjacency list)
    # =========================================================================
    try:
        os.makedirs(args.outdir)
    except FileExistsError:
        raise Exception('{} directory already exists.'.format(args.outdir) +
                        'Please specify another location or rename it.')

    N1 = float_matrix.shape[0]
    matrix = binarize_matrix(args, float_matrix)

    # ++++ Tracking clust/un-clustered bits to avoid re-computations ++++++++++
    unclust_bit = ba(N1)
    unclust_bit.setall(1)
    clustered_bit = unclust_bit.copy()
    clustered_bit.setall(0)
    zeros = np.zeros(N1, dtype=np.int32)

    # ++++ Save clusters in an array (1 .. N) +++++++++++++++++++++++++++++++++
    clusters_array = np.zeros(N1, dtype=np.int32)
    NCLUSTER = 0
    clustered = set()
    n_members = []

    # ++++ Coloring ordered vertices (1 .. N) +++++++++++++++++++++++++++++++++
    degrees = calc_matrix_degrees(unclust_bit, matrix)
    ordered_by_degs = degrees.argsort()[::-1]
    colors = colour_matrix(ordered_by_degs, matrix)

    # =========================================================================
    # 2. Main algorithm: BitQT !
    # =========================================================================
    while any(degrees):
        NCLUSTER += 1
        # ++++ Find a big clique early ++++++++++++++++++++++++++++++++++++++++
        big_node = degrees.argmax()
        bit_clique, big_clique = do_bit_cascade(big_node, degrees, colors,
                                                matrix, 0)
        big_clique_size = big_clique.size
        # ++++ Find promising nodes +++++++++++++++++++++++++++++++++++++++++++
        biggers = degrees > big_clique_size
        biggers[big_clique] = False
        cluster_colors = colors[big_clique]
        biggers_colors = colors[biggers]
        promising_colors = np.setdiff1d(biggers_colors, cluster_colors)
        promising_nodes = deque()
        for x in promising_colors:
            promising_nodes.extend(((colors == x) & biggers).nonzero()[0])
        # ++++ Explore all promising nodes ++++++++++++++++++++++++++++++++++++
        cum_found = big_clique
        while promising_nodes:
            node = promising_nodes.popleft()
            try:
                bit_clique, clique = do_bit_cascade(node, degrees, colors,
                                                    matrix, big_clique_size)
                CLIQUE_SIZE = len(clique)
            except TypeError:
                CLIQUE_SIZE = 0
            # ++++ Cumulative update only if biggers candidates are found +++++
            if CLIQUE_SIZE > big_clique_size:
                big_node = node
                big_clique = clique
                big_clique_size = big_clique.size
                # ++++ Repeat previous condition ++++++++++++++++++++++++++++++
                cum_found = np.concatenate((cum_found, big_clique))
                biggers = degrees > big_clique_size
                biggers[cum_found] = False
                cluster_colors = colors[big_clique]
                biggers_colors = colors[biggers]
                promising_colors = np.setdiff1d(biggers_colors, cluster_colors)
                promising_nodes = deque()
                for x in promising_colors:
                    promising_nodes.extend(
                        ((colors == x) & biggers).nonzero()[0])
        n_members.append(big_clique_size)
        if (big_clique_size < args.min_clust_size) or (
                NCLUSTER == args.clust_id):
            break

        # ++++ Save new cluster & update NCLUSTER +++++++++++++++++++++++++++++
        clusters_array[big_clique] = NCLUSTER
        # ++++ Update (un)clustered_bit +++++++++++++++++++++++++++++++++++++++
        clustered.update(big_clique)
        clustered_bit = set_to_bitarray(clustered, N1)
        unclust_bit = ~clustered_bit
        # ++++ Hard erasing of clustered frames from matrix +++++++++++++++++++
        degrees = zeros.copy()
        for x in unclust_bit[:N1].itersearch(ba('1')):
            degrees[x] = matrix[x].count()
            if bu.count_and(matrix[x], clustered_bit):
                matrix[x] &= (matrix[x] ^ clustered_bit)

    # =========================================================================
    # 3. Output
    # =========================================================================
    # saving pickle for api debugging tests
    # outname = os.path.basename(args.topology).split('.')[0]
    # pickle_to_file(clusters_array, os.path.join(args.outdir,
    #                                             '{}.pick'.format(outname)))
    # saving VMD visualization script
    # to_VMD(args.outdir, args.topology, args.first, args.last, N1, args.stride,
    #        clusters_array[:N1])
    # saving clustering info  files
    frames_stats = get_frames_stats(N1, args.first, args.last, args.stride,
                                    clusters_array[:N1], args.outdir)
    cluster_stats = get_cluster_stats(clusters_array[:N1], args.outdir)
    print('\n\nNormal Termination of BitQT :)')
    return frames_stats, cluster_stats


# =============================================================================
#
# =============================================================================
args = argparse.Namespace()
args.first = 0
args.stride = 1
args.min_clust_size = 2
args.cutoff = 3
args.n_clust = np.inf
args.outdir = 'bitQT_outputs'

sampling_path = '/home/roy.gonzalez-aleman/RoyHub/stdock/data/autodock/results_1x2000/docking_poses/poses_coords.pdb'
sampling = prd.parsePDB(sampling_path)
ensemble = prd.Ensemble()
ensemble.setAtoms(sampling)
ensemble.setCoords(sampling.getCoords())
[ensemble.addCoordset(cs) for cs in sampling.getCoordsets()]
rmsd_matrix = ensemble.getRMSDs(pairwise=True)
args.last = len(rmsd_matrix)

frames, clusts = bitqt_matrix(args, rmsd_matrix)

import matplotlib.pyplot as plt
dejavu = set()
data = []
for frame in range(frames.shape[0]):
    clust_id = frames.iloc[frame].cluster_id
    if clust_id not in dejavu:
        data.append((frame, clust_id))
        dejavu.add(clust_id)

simple_coverage_x = []
simple_coverage_y = []
n = 0
for d in data:
    simple_coverage_x.append(d[0])
    simple_coverage_y.append(n + 1)
    n += 1

plt.plot(simple_coverage_x, simple_coverage_y)
plt.show()
