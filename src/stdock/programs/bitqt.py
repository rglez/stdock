# Created by roy.gonzalez-aleman at 03/12/2023
import os
import pickle
from collections import deque, OrderedDict

import mdtraj as md
import numpy as np
import pandas as pd
from bitarray import bitarray as ba
from bitarray import util as bu

valid_tops = {'pdb', 'pdb.gz', 'h5', 'lh5', 'prmtop', 'parm7', 'prm7', 'psf',
              'mol2', 'hoomdxml', 'gro', 'arc', 'hdf5', 'gsd'}

valid_trajs = {'arc', 'dcd', 'binpos', 'xtc', 'trr', 'hdf5', 'h5', 'ncdf',
               'netcdf', 'nc', 'pdb.gz', 'pdb', 'lh5', 'crd', 'mdcrd',
               'inpcrd', 'restrt', 'rst7', 'ncrst', 'lammpstrj', 'dtr', 'stk',
               'gro', 'xyz.gz', 'xyz', 'tng', 'xml', 'mol2', 'hoomdxml', 'gsd'}


def binarize_matrix(float_matrix, cutoff):
    N = float_matrix.shape[0]
    cutoff = np.full(N, cutoff, dtype=np.float32)
    matrix = dict()
    for i in range(N):
        rmsd_ = float_matrix[i]
        vector_np = np.less_equal(rmsd_, cutoff)
        bit_array = ba()
        bit_array.pack(vector_np.tobytes())
        bit_array.fill()
        matrix.update({i: bit_array})
    return matrix


def is_valid_traj(traj, valid_trajs):
    traj_ext = traj.split('.')[-1]
    if traj_ext not in valid_trajs:
        raise ValueError('The trajectory format "{}" '.format(traj_ext) +
                         'is not available. Valid trajectory formats '
                         'are: {}'.format(valid_trajs))
    return True


def traj_needs_top(traj):
    traj_ext = traj.split('.')[-1]
    if traj_ext in ['h5', 'lh5', 'pdb']:
        return False
    return True


def is_valid_top(topology, valid_tops):
    try:
        top_ext = topology.split('.')[-1]
    except AttributeError:
        raise ValueError('You should pass a topology object. '
                         'Valid topology formats are: {}'.format(valid_tops))

    if top_ext not in valid_tops:
        raise ValueError('The topology format "{}"'.format(top_ext) +
                         'is not available. Valid topology formats'
                         'are: {}'.format(valid_tops))
    return True


def load_raw_traj(traj, valid_trajs, topology=None):
    if is_valid_traj(traj, valid_trajs) and traj_needs_top(traj):
        if is_valid_top(topology, valid_tops):
            return md.load(traj, top=topology)

    if is_valid_traj(traj, valid_trajs) and not traj_needs_top(traj):
        return md.load(traj)


def shrink_traj_selection(traj, selection):
    if selection != 'all':
        try:
            sel_indx = traj.topology.select(selection)
        except Exception:
            raise ValueError('Specified selection is invalid')
        if sel_indx.size == 0:
            raise ValueError('Specified selection corresponds to no atoms')
        traj = traj.atom_slice(sel_indx, inplace=True)
    return traj


def shrink_traj_range(first, last, stride, traj):
    # Calculate range of available intervals ----------------------------------
    N = traj.n_frames
    first_range = range(0, N - 1)
    last_range = range(first + 1, N)
    try:
        delta = last - first
    except TypeError:
        delta = N - first
    stride_range = range(1, delta)
    # Raising if violations ---------------------------------------------------
    if first not in first_range:
        raise ValueError('"first" parameter should be in the interval [{},{}]'
                         .format(first_range.start, first_range.stop))
    if last and (last not in last_range):
        raise ValueError('"last" parameter should be in the interval [{},{}]'
                         .format(last_range.start, last_range.stop))
    if stride not in stride_range:
        raise ValueError('"stride" parameter should be in the interval [{},{}]'
                         .format(stride_range.start, stride_range.stop))
    # Slicing trajectory ------------------------------------------------------
    sliced = slice(first, last, stride)
    if sliced not in [slice(0, N, 1), slice(0, None, 1)]:
        return traj.slice(sliced)
    return traj


def calc_rmsd_matrix(trajectory, cutoff):
    trajectory.center_coordinates()
    cutoff = np.full(trajectory.n_frames, cutoff / 10, dtype=np.float32)
    matrix = OrderedDict()
    to_explore = range(trajectory.n_frames)
    for i in to_explore:
        rmsd_ = md.rmsd(trajectory, trajectory, i, precentered=True)
        vector_np = np.less_equal(rmsd_, cutoff)
        bit_array = ba()
        bit_array.pack(vector_np.tobytes())
        bit_array.fill()
        matrix.update({i: bit_array})
    return matrix


def calc_matrix_degrees(unclustered_bit, matrix):
    one = ba('1')
    degrees = np.zeros(len(unclustered_bit), dtype=np.int32)
    for node in unclustered_bit.itersearch(one):
        try:
            degrees[node] = matrix[node].count()
        except KeyError:
            pass
    return degrees


def colour_matrix(degrees, matrix):
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
        colored.pack(colors.astype(bool).tobytes())
        if colored.count(0) == 0:
            break
    return colors


def bitarray_to_np(bitarr):
    return np.unpackbits(bitarr).astype(bool)


def do_bit_cascade(big_node, degrees, colors, matrix, max_):
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
    zero_arr = np.zeros(N, dtype=bool)
    zero_arr[list(set_)] = 1
    bitarr = ba()
    bitarr.pack(zero_arr.tobytes())
    return bitarr


def pickle_to_file(data, file_name):
    with open(file_name, 'wb') as file:
        pickle.dump(data, file)
    return file_name


def top_has_coords(topology):
    try:
        tt = md.load(topology)
    except OSError:
        return False
    return tt.xyz.shape[0]


def to_VMD(outdir, topology, first, N1, last, stride, final_array):
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


def get_frames_stats(N1, clusters, outdir):
    frames_df = pd.DataFrame(columns=['frame', 'cluster_id'])
    frames_df['frame'] = range(N1)
    frames_df['cluster_id'] = clusters
    with open(os.path.join(outdir, 'frames_statistics.txt'), 'wt') as on:
        frames_df.to_string(buf=on, index=False)
    return frames_df


def get_cluster_stats(clusters, outdir):
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


def bitqt_matrix(matrix, min_clust_size, max_num_clust, outdir):
    # =========================================================================
    # 1. Creating binary matrix (adjacency list)
    # =========================================================================
    try:
        os.makedirs(outdir)
    except FileExistsError:
        raise Exception('{} directory already exists.'.format(outdir) +
                        'Please specify another location or rename it.')

    N1 = len(matrix)
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
        if (big_clique_size < min_clust_size) or (
                NCLUSTER == max_num_clust):
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
    # outname = os.path.basename(topology).split('.')[0]
    # pickle_to_file(clusters_array, os.path.join(outdir,
    #                                             '{}.pick'.format(outname)))
    # saving VMD visualization script
    # to_VMD(outdir, topology, first, last, N1, stride,
    #        clusters_array[:N1])
    # saving clustering info  files
    frames_stats = get_frames_stats(N1, clusters_array[:N1], outdir)
    cluster_stats = get_cluster_stats(clusters_array[:N1], outdir)
    print('\n\nNormal Termination of BitQT :)')
    return frames_stats, cluster_stats


def get_matrix_from_topotraj(trajectory, topology, selection, cutoff):
    trajectory = load_raw_traj(trajectory, valid_trajs, topology)
    trajectory = shrink_traj_selection(trajectory, selection)
    trajectory.center_coordinates()
    matrix = calc_rmsd_matrix(trajectory, cutoff)
    return matrix




# =============================================================================
# debugging area
# =============================================================================
import prody as prd

trajectory = '/home/roy.gonzalez-aleman/RoyHub/stdock/scripts/dude_selection/01_conformers_generation/nos1_conformers.pdb'
topology = '/home/roy.gonzalez-aleman/RoyHub/stdock/scripts/dude_selection/01_conformers_generation/top_nos1_conformers.pdb'
selection = 'all'
cutoff = 2

ensemble = prd.Ensemble(prd.parsePDB(trajectory))
float_matrix = ensemble.getRMSDs(pairwise=True)

bin_matrix_topotraj = get_matrix_from_topotraj(trajectory, topology, selection, cutoff)
bin_matrix_floatmatrix = binarize_matrix(float_matrix, cutoff)

bitqt_matrix(bin_matrix_topotraj, 2, 1000000, './topotraj')

bitqt_matrix(bin_matrix_floatmatrix, 2, 1000000, './floatmatrix')
