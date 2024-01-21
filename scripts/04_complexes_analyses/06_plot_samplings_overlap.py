# Created by roy.gonzalez-aleman at 20/01/2024
import commons as cmn

# =============================================================================
# User-defined parameters
# =============================================================================
proj_dir = cmn.proj_dir
top_dir = '/scripts/04_complexes_analyses/03_clusterings'
mappings_dir = '/scripts/04_complexes_analyses/01_rec_coverage'
out_dir = '/scripts/04_complexes_analyses/04_plot_clusterings'
os.makedirs(out_dir, exist_ok=True)

top_x = 20
cases = [join(top_dir, case) for case in os.listdir(top_dir)]
data_dict = cmn.recursive_defaultdict()

# ==== Plot
for case_top_dir in cases:
    case = basename(case_top_dir)
    clusters_no_super = cmn.unpickle_from_file(
        next(cmn.recursive_finder('*_nosuper*.pick', case_top_dir)))
    clusters_super = cmn.unpickle_from_file(
        next(cmn.recursive_finder('*_super*.pick', case_top_dir)))
    labels_ids = np.load(
        next(cmn.recursive_finder('*_identifiers*', case_top_dir)))
