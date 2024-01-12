# Created by roy.gonzalez-aleman at 09/01/2024
import itertools as it
import os
from collections import Counter, deque
from os.path import join

import matplotlib.transforms as transforms
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

import commons as cmn

# todo: put the clusters around the asr in top of the figure?
# =============================================================================
# User-defined parameters
# =============================================================================
# top_dir = '/home/roy.gonzalez-aleman/RoyHub/stdock/scripts/dude_analyses/03_clusterings'
# mappings_dir = '/home/roy.gonzalez-aleman/RoyHub/stdock/scripts/dude_analyses/01_rec_coverage'
# top_x = 20
#
# cases = [join(top_dir, case) for case in os.listdir(top_dir)]
# data_dict = cmn.recursive_defaultdict()
# for case_top_dir in cases:
#     clusters_no_super = cmn.unpickle_from_file(
#         next(cmn.recursive_finder('*_nosuper*.pick', case_top_dir)))
#     clusters_super = cmn.unpickle_from_file(
#         next(cmn.recursive_finder('*_super*.pick', case_top_dir)))
#     labels_ids = np.load(
#         next(cmn.recursive_finder('*_identifiers*', case_top_dir)))
#
#     # =========================================================================
#     # Clusters no_super (internal sampling)
#     # =========================================================================
#     # Setting the kind of clustering to deal with
#     clusters = clusters_no_super

    # Creating & sorting a df with cluster counts
    labels_only = np.asarray(['_'.join(x.split('_')[:-1]) for x in labels_ids])
    counters = [Counter(labels_only[x]) for x in clusters]
    df_counts = pd.DataFrame(counters)

    mapping = next(cmn.recursive_finder('*mapping*', mappings_dir))
    columns = pd.read_table(mapping, sep='\s+').columns
    df_counts = df_counts[columns].T

    # Restrict to a top-x fraction of clusters
    ranks_only = np.asarray([int(x.split('_')[-1]) for x in labels_ids])
    seeds_ranks = [ranks_only[x[0]] for x in clusters]
    seeds_progs = [labels_only[x[0]] for x in clusters]
    top_seeds_ranks = list(it.takewhile(lambda x: x <= top_x, seeds_ranks))
    N = len(top_seeds_ranks)
    top_x_counts = df_counts[range(N)]
    top_seeds_prog = seeds_progs[:N]

    # Request info for plotting
    clusters_tot_poses = df_counts.sum()
    programs_tot_poses = top_x_counts.sum(axis=1)
    percent = df_counts.divide(df_counts.sum(axis=1), axis=0) * 100

    best_top = [int(labels_ids[cluster[0]].split('_')[-1]) for cluster in
                clusters]
    summation_right = (top_x_counts > 0).sum(axis=1) / N * 100
    M = top_x_counts.shape[0]
    summation_top = (top_x_counts > 0).sum() / M * 100

    # =========================================================================
    #
    # =========================================================================
    # ==== General formatting
    cmn.reset_matplotlib()
    cmn.generic_matplotlib(width=(3.25, 5))
    cmap_original = plt.get_cmap('Greys')
    cmap = cmn.chop_cmap_frac(cmap_original, 0.2)
    cmap.set_under(color='white')
    color_1 = 'gray'
    size_1 = 8
    alpha_1 = 0.75

    # ==== Layout specifications
    fig = plt.figure(layout="constrained")
    fig.tight_layout()
    gs = GridSpec(13, 13, figure=fig, wspace=0.25, hspace=0.)
    central = fig.add_subplot(gs[1:, :-3])
    # top = fig.add_subplot(gs[:1, :-3], sharex=central)
    right = fig.add_subplot(gs[1:, -3:-1], sharey=central)
    color_bar = fig.add_subplot(gs[1:, -1:])

    # ==== Plotting
    # 1. Central Data
    mappings = central.pcolor(top_x_counts, cmap=cmap, vmin=0.5)

    # 2. Right Coverage
    positions_right = [x + 0.5 for x in range(len(summation_right))]
    right.barh(positions_right, summation_right, color=color_1,
               alpha=0.5, lw=0.2, height=1, edgecolor=color_1)

    # 3. Color Bar
    cbar = fig.colorbar(mappings, ax=central, cax=color_bar, cmap=cmap)

    # 4. Top ASR

    # ==== Labels
    central.set_xlabel('Top-ranked clusters', weight='bold')
    right.set_xlabel('% Coverage', weight='bold')
    cbar.set_label('# members', weight='bold')

    # ==== Ticks

    # 2. Label indexing
    final_names = cmn.get_final_names(top_x_counts.index)

    # 3. Putting text
    df = pd.DataFrame(final_names, columns=['prog', 'sf', 'exh'])
    bicolor = deque(['k', color_1])
    trans = transforms.blended_transform_factory(central.transAxes,
                                                 central.transData)

    prog_indices = df.groupby(['prog']).indices
    prog_colors = {}
    for label in prog_indices:
        color = bicolor[0]
        [prog_colors.update({x: color}) for x in prog_indices[label]]
        bicolor.rotate()
        pos = prog_indices[label].mean()
        central.text(-0.35, pos + 0.5, str(label),
                     transform=trans,
                     horizontalalignment='center',
                     verticalalignment='center', rotation=90,
                     fontdict={'weight': 'bold', 'size': 6,
                               'color': prog_colors[int(pos)]})

    sf_indices = df.groupby(['prog', 'sf']).indices
    for label in sf_indices:
        sf_name = label[-1]
        pos = sf_indices[label].mean()
        central.text(-0.2, pos + 0.5, sf_name,
                     transform=trans,
                     horizontalalignment='center',
                     verticalalignment='center', rotation=0,
                     fontdict={'weight': 'normal', 'size': 6,
                               'color': prog_colors[int(pos)]})

    exh_indices = df.groupby('exh').indices
    for label in exh_indices:
        for pos in exh_indices[label]:
            central.text(-0.05, pos + 0.5, str(label),
                         transform=trans,
                         horizontalalignment='center',
                         verticalalignment='center', rotation=0,
                         fontdict={'weight': 'normal', 'size': 6,
                                   'color': prog_colors[int(pos)]})

    # ==== Grids
    [central.axhline(x, ls='--', lw=0.5, c=color_1, alpha=alpha_1)
     for x in range(M)]
    [central.axhline(x, ls='-', lw=0.5, c='k')
     for x in range(0, M, 3)]
    [right.axhline(x, ls='-', lw=0.5, c=color_1)
     for x in range(0, M, 3)]
    [right.axvline(x, color='gray', ls='--', lw=0.5) for x in [25, 50, 75]]
    central.grid(axis='x', ls='--', color='k', lw=0.5)

    # ==== Hidings & Reorderings
    # plt.setp(top.get_yticklabels(), visible=False)
    # plt.setp(top.get_xticklabels(), visible=False)
    # top.spines['right'].set_visible(False)
    # top.spines['left'].set_visible(False)
    # top.xaxis.set_label_position('top')
    # top.set_yticks([])
    # top.spines['top'].set_visible(False)

    plt.setp(right.get_yticklabels(), visible=False)
    right.set_yticks([])
    plt.setp(right.get_yticklabels(), visible=False)
    right.tick_params(axis='y', which='both', left=False)
    right.spines['top'].set_visible(True)
    right.xaxis.set_label_position('top')
    right.xaxis.set_ticks_position('top')

    # ==== Saving
    plt.savefig(join(case_top_dir, f'clusters_no_super.png'),
                bbox_inches='tight')
    plt.close()
