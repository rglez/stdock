# Created by roy.gonzalez-aleman at 10/03/2024
import os
from collections import deque
from os.path import dirname, join, basename

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import ticker
import commons as cmn
import matplotlib.colors as colors

# ==== User-specified arguments ===============================================
insights_pickle = ('/home/roy.gonzalez-aleman/RoyHub/stdock/scripts/'
                   '04_complexes_analyses/07_insights/insights.pickle')
num_rot_file = 'scripts/01_complexes_selection/03_selection/report_selection.txt'
out_dir = dirname(insights_pickle)


# =============================================================================

def get_data(string, scoring_functions, cases_rot):
    data = cmn.recursive_defaultdict()
    for sf in scoring_functions:
        for case in cases_rot:
            try:
                data[sf][case] = all_info[case][sf][string]
            except KeyError:
                pass
    return pd.DataFrame(data).T


insights = cmn.unpickle_from_file(insights_pickle)

all_info = {}
for case in insights.keys():
    out_dir = join(out_dir, case)
    os.makedirs(out_dir, exist_ok=True)

    # Individual info update
    info = dict()
    for program in ['plants', 'autodock4', 'vina', 'smina', 'gnina', 'qvinaw']:
        for sf in sorted(insights[case][program].keys()):
            for level in ['low', 'medium', 'high']:
                mapping_name = f'{program}_{sf}_{level}'
                info.update({mapping_name: {}})

                n_poses = insights[case][program][sf][level]['n_poses']
                coverage = insights[case][program][sf][level]['coverage']
                clust_int = insights[case][program][sf][level]['n_clust_int']
                clust_ext = insights[case][program][sf][level]['n_clust_ext']
                time = insights[case][program][sf][level]['time']
                ram = insights[case][program][sf][level]['ram']

                info[mapping_name]['n_poses'] = n_poses
                info[mapping_name]['coverage'] = coverage
                info[mapping_name]['n_clust_int'] = clust_int
                info[mapping_name]['n_clust_ext'] = clust_ext
                info[mapping_name]['time'] = time
                info[mapping_name]['ram'] = ram
                all_info.update({case: info})

# ==== Cases
names = ['n_rot', 'path']
raw = pd.read_table(num_rot_file, sep=':', header=None, names=names)
basename_rot = dict(zip([basename(x) for x in raw.path], raw.n_rot))
s_functions = list(info.keys())
cases_rot = {x.split('_')[1]: basename_rot[x] for x in basename_rot}
cases_labels = [f'{x} ({cases_rot[x]: >2})' for x in cases_rot if
                x in all_info]

# %%==== Data to plot
data_plot = {0: get_data('n_poses', s_functions, cases_rot),
             3: get_data('coverage', s_functions, cases_rot),
             1: get_data('n_clust_int', s_functions, cases_rot),
             4: get_data('n_clust_ext', s_functions, cases_rot),
             2: get_data('time', s_functions, cases_rot) / 60,
             5: get_data('ram', s_functions, cases_rot)}

# ==== General formatting
n_cases = len(cases_rot)
n_sf = len(s_functions)

cmn.reset_matplotlib()
cmn.generic_matplotlib(width=(14, 9))

color_1 = 'grey'
color_2 = 'grey'
alpha_1 = 1
cmap = 'viridis'

# ==== Layout specifications
fig = plt.figure()
# fig.tight_layout()
gs = GridSpec(2, 3, figure=fig, wspace=0.08, hspace=0.14)

zero = fig.add_subplot(gs[:1, :1])
one = fig.add_subplot(gs[:1, 1:2], sharey=zero)
two = fig.add_subplot(gs[:1, 2:], sharey=zero)

three = fig.add_subplot(gs[1:, :1])
four = fig.add_subplot(gs[1:, 1:2], sharey=three)
five = fig.add_subplot(gs[1:, 2:], sharey=three)

# ==== Plotting
poses = zero.pcolormesh(data_plot[0], cmap=cmap, edgecolors=color_1, lw=0.1,
                        norm=colors.LogNorm())
clustin = one.pcolormesh(data_plot[1], cmap=cmap, edgecolors=color_1, lw=0.1,
                         norm=colors.LogNorm())
timeit = two.pcolormesh(data_plot[2], cmap=cmap, edgecolors=color_1, lw=0.1,
                        norm=colors.LogNorm())
covered = three.pcolormesh(data_plot[3], cmap=cmap, edgecolors=color_1, lw=0.1,
                           norm=colors.LogNorm())
clustext = four.pcolormesh(data_plot[4], cmap=cmap, edgecolors=color_1, lw=0.1,
                           norm=colors.LogNorm())
ramit = five.pcolormesh(data_plot[5], cmap=cmap, edgecolors=color_1, lw=0.1,
                        norm=colors.LogNorm())

# ==== Labels
four.set_xlabel('Benchmarked Complexes ID (# rot bonds)', weight='bold',
                size=16, labelpad=12)
zero.set_title('# Lig Poses', fontsize=15)
three.set_title('Receptor Coverage [%]', fontsize=15)

one.set_title('# Clusters (super)', fontsize=15)
four.set_title('# Clusters (no-super)',  fontsize=15)

two.set_title('Runtime [min]', fontsize=15)
five.set_title('RAM [GB]', fontsize=15)

# 2. Label indexing
col_names = list(data_plot[0].index)
final_names = cmn.get_final_names(col_names)

# ==== Ticks

# 3. Putting text
df = pd.DataFrame(final_names, columns=['prog', 'sf', 'exh'])
bicolor = deque(['k', 'darkgrey'])

prog_indices = df.groupby(['prog']).indices

for selected in [zero, three]:
    prog_colors = {}
    for label in ['PLANTS', 'AD4', 'VINA', 'SMINA', 'GNINA', 'QVINAW']:
        color = bicolor[0]
        [prog_colors.update({x: color}) for x in prog_indices[label]]
        bicolor.rotate()
        pos = prog_indices[label].mean()
        selected.text(-7.5, pos + 0.5, str(label),
                      horizontalalignment='center',
                      verticalalignment='center', rotation=0,
                      fontdict={'weight': 'bold', 'size': 10,
                                'color': prog_colors[int(pos)]})

    sf_indices = df.groupby(['prog', 'sf']).indices
    for label in sf_indices:
        sf_name = label[-1]
        pos = sf_indices[label].mean()
        selected.text(-3.5, pos + 0.5, sf_name,
                      horizontalalignment='center',
                      verticalalignment='center', rotation=0,
                      fontdict={'weight': 'normal', 'size': 10,
                                'style': 'italic',
                                'color': prog_colors[int(pos)]})

    exh_indices = df.groupby('exh').indices
    for label in exh_indices:
        for pos in exh_indices[label]:
            selected.text(-1, pos + 0.5, str(label),
                          horizontalalignment='center',
                          verticalalignment='center', rotation=90,
                          fontdict={'weight': 'normal', 'size': 8,
                                    'color': prog_colors[int(pos)]})

displaced = [x + 0.5 for x in range(len(cases_labels))]

# ==== Grids
M = len(col_names)
for i, axis in enumerate([zero, one, two, three, four, five]):
    [axis.axhline(x, ls='-', lw=1, c=color_1)
     for x in range(0, M, 3)]

# ==== Hidings & Reorderings
zero.set_yticks([])

three.set_yticks([])
for x in [three, four, five]:
    x.set_xticks(displaced, cases_labels, rotation=90, size=9)
for x in [zero, one, two]:
    x.set_xticks([])

# ==== Colorbars
cbar0 = fig.colorbar(poses, ax=zero, cmap=cmap, pad=0.01)
cbar1 = fig.colorbar(clustin, ax=one, cmap=cmap, pad=0.01)
cbar2 = fig.colorbar(timeit, ax=two, cmap=cmap, pad=0.01)

# for cbar in [cbar0, cbar1, cbar2]:
#     cbar.locator = ticker.MaxNLocator(nbins=5)
#     cbar.update_ticks()
#     cbar.formatter.set_powerlimits((0, 0))
#     cbar.formatter.set_useMathText(True)

cbar = fig.colorbar(covered, ax=three, cmap=cmap, pad=0.01)
cbar = fig.colorbar(clustext, ax=four, cmap=cmap, pad=0.01)
cbar = fig.colorbar(ramit, ax=five, cmap=cmap, pad=0.01)

# plt.show()
plt.savefig('benchmark_test.png', bbox_inches='tight')
plt.close()
