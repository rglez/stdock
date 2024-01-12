# Created by roy.gonzalez-aleman at 02/01/2024
from collections import deque
from os.path import join, split

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec

import commons as cmn

# =============================================================================
# User-defined parameters
# =============================================================================
top_dir = '/home/roy.gonzalez-aleman/RoyHub/stdock/scripts/dude_analyses/01_rec_coverage'
mappings = list(cmn.recursive_finder('mappings.txt', top_dir))

for mapping in mappings:
    out_dir = split(mapping)[0]

    # ==== Getting data
    mapping_df = pd.read_table(mapping, sep='\s+')
    total_res_nums, total_num_sf = mapping_df.shape
    summation_top = (mapping_df.T > 0).sum() / total_num_sf * 100
    summation_right = (mapping_df > 0).sum() / total_res_nums * 100
    col_names = list(mapping_df.T.index)
    M = len(col_names)
    # todo: get the ASR as  within 5 of resname LIG from PDB
    asr = [185, 115, 170, 117, 171, 182, 111, 168, 114, 186, 181, 148, 116,
           113, 149, 169, 164, 165, 183, 184, 166, 112, 167, 180]

    # ==== General formatting
    cmn.reset_matplotlib()
    cmn.generic_matplotlib(width=(3.25, 5))
    cmap_original = plt.get_cmap('Greys')
    cmap = cmn.chop_cmap_frac(cmap_original, 0.1)
    cmap.set_under(color='white')
    color_1 = 'gray'
    color_2 = 'gray'
    size_1 = 8
    alpha_1 = 0.75

    # ==== Layout specifications
    fig = plt.figure(layout="constrained")
    fig.tight_layout()
    gs = GridSpec(13, 13, figure=fig, wspace=0.25, hspace=0.)
    central = fig.add_subplot(gs[1:, :-3])
    top = fig.add_subplot(gs[:1, :-3], sharex=central)
    right = fig.add_subplot(gs[1:, -3:-1], sharey=central)
    color_bar = fig.add_subplot(gs[1:, -1:])

    # ==== Plotting
    # 1. Central Data
    mappings = central.pcolor(mapping_df.T, cmap=cmap, vmin=0.5)

    # 2. Right Coverage
    positions_right = [x + 0.5 for x in range(len(summation_right))]
    right.barh(positions_right, summation_right, color=color_1,
               alpha=0.5, lw=0.2, height=1, edgecolor=color_1)

    # 3. Color Bar
    cbar = fig.colorbar(mappings, ax=central, cax=color_bar, cmap=cmap)

    # 4. Top ASR
    [top.axvline(x, ymin=0., ymax=0.5, color=color_1) for x in asr]

    # ==== Labels
    central.set_xlabel('Protein residue number', weight='bold')
    right.set_xlabel('% Coverage', weight='bold')
    cbar.set_label('# close ligand poses', weight='bold')

    # ==== Ticks

    # 2. Label indexing
    final_names = cmn.get_final_names(col_names)

    # 3. Putting text
    df = pd.DataFrame(final_names, columns=['prog', 'sf', 'exh'])
    bicolor = deque(['k', color_1])

    prog_indices = df.groupby(['prog']).indices
    prog_colors = {}
    for label in prog_indices:
        color = bicolor[0]
        [prog_colors.update({x: color}) for x in prog_indices[label]]
        bicolor.rotate()
        pos = prog_indices[label].mean()
        central.text(-150, pos + 0.5, str(label),
                     horizontalalignment='center',
                     verticalalignment='center', rotation=90,
                     fontdict={'weight': 'bold', 'size': 6,
                               'color': prog_colors[int(pos)]})

    sf_indices = df.groupby(['prog', 'sf']).indices
    for label in sf_indices:
        sf_name = label[-1]
        pos = sf_indices[label].mean()
        central.text(-88, pos + 0.5, sf_name,
                     horizontalalignment='center',
                     verticalalignment='center', rotation=0,
                     fontdict={'weight': 'normal', 'size': 6,
                               'style': 'italic',
                               'color': prog_colors[int(pos)]})

    exh_indices = df.groupby('exh').indices
    for label in exh_indices:
        for pos in exh_indices[label]:
            central.text(-20, pos + 0.5, str(label),
                         horizontalalignment='center',
                         verticalalignment='center', rotation=0,
                         fontdict={'weight': 'normal', 'size': 6,
                                   'color': prog_colors[int(pos)]})

    # ==== Grids
    [central.axhline(x, ls='--', lw=0.5, c=color_2, alpha=alpha_1)
     for x in range(M)]
    [central.axhline(x, ls='-', lw=0.5, c=color_2)
     for x in range(0, M, 3)]
    [right.axhline(x, ls='-', lw=0.5, c=color_1)
     for x in range(0, M, 3)]
    [right.axvline(x, color='gray', ls='--', lw=0.5) for x in [25, 50, 75]]
    central.grid(axis='x', ls='--', color='k', lw=0.5)

    # ==== Hidings & Reorderings
    plt.setp(top.get_yticklabels(), visible=False)
    plt.setp(top.get_xticklabels(), visible=False)
    top.spines['right'].set_visible(False)
    top.spines['left'].set_visible(False)
    top.xaxis.set_label_position('top')
    top.set_yticks([])
    top.spines['top'].set_visible(False)

    plt.setp(right.get_yticklabels(), visible=False)
    right.set_yticks([])
    plt.setp(right.get_yticklabels(), visible=False)
    right.tick_params(axis='y', which='both', left=False)
    right.spines['top'].set_visible(True)
    right.xaxis.set_label_position('top')
    right.xaxis.set_ticks_position('top')

    # ==== Saving
    plt.savefig(join(out_dir, f'rec_surface_coverage.png'),
                bbox_inches='tight')
    plt.close()
