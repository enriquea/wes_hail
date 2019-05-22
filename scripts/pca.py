"""
Run population PCA on QCed and AF-filtered MatrixTable.
Perform variant LD-pruning before run PCA.

"""

import argparse
import time
from config import HAIL_LOG_PATH
from utils.hail_functions import *


def main(args):
    # Initializing Hail on cluster mode
    init_hail_on_cluster(app_name='pca-workflow',
                         tmp_dir='/mnt/nfs/mdatanode/hail-temp',
                         log_file=HAIL_LOG_PATH,
                         local_mode=False)

    # Read MatrixTable
    mt = hl.read_matrix_table(args.mt_path)

    # Filter by AF (keep common variants)
    mt_common = (mt
                 .filter_rows(mt.variant_qc.AF[1] >= args.af_threshold, keep=True)
                 )

    if args.ld_pruning:
        # LD pruning
        # Avoid filtering / missingness entries (genotypes) before run LP pruning
        # Zulip Hail support issue -> "BlockMatrix trouble when running pc_relate"
        mt_common = mt_common.unfilter_entries()

        pruned_variant_table = hl.ld_prune(mt_common.GT,
                                           r2=args.r2,
                                           bp_window_size=args.bp_window_size,
                                           memory_per_core=512)

        # Keep LD-pruned variants
        mt_common = (mt_common
                     .filter_rows(hl.is_defined(pruned_variant_table[mt_common.locus, mt_common.alleles]), keep=True)
                     )

    # Run PCA using hwe normalized method from Hail.
    if args.export_loadings:
        eigenvalues, pc_scores, loadings = hl.hwe_normalized_pca(mt_common.GT, k=args.n_pcs, compute_loadings=True)
    else:
        eigenvalues, pc_scores, _ = hl.hwe_normalized_pca(mt_common.GT, k=args.n_pcs)

    # getting number of samples/variants
    n_samples, n_variants = mt_common.count()

    # Define global annotations
    date = time.strftime("%d%m%Y")
    annotations = {'date': date,
                   'input_path': args.mt_path,
                   'af_cutoff': args.af_threshold,
                   'eigenvalues': eigenvalues,
                   'nummber_of_pcs': args.n_pcs,
                   'pca_method': 'hail.hwe_normalized',
                   'n_samples': n_samples,
                   'n_variants': n_variants,
                   'ld_pruned': args.ld_pruning,
                   'r2': args.r2}

    ht = add_global_annotations(ds=pc_scores, annotations=annotations, overwrite=True)

    # Export PCs table
    (ht
     .write(output=f'{args.ht_output_path}_scores_{date}.ht')
     )

    # Export loadings
    if args.export_loadings:
        loadings.write(output=f'{args.ht_output_path}_loadings_{date}.ht')

    # Write PCs table to file
    if args.write_to_file:
        # Unpacking nested PCs table scores
        pca_table = (pc_scores
                     .annotate(**{'PC' + str(k + 1): pc_scores.scores[k] for k in range(0, args.n_pcs)})
                     .drop('scores')
                     )

        # Export table to file
        pca_table.export(output=f'{args.ht_output_path}pc_scores_{date}.tsv')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--mt_path', help='Path to MatrixTable with computed QC metrics', type=str, default=None)
    parser.add_argument('--ht_output_path', help='Output HailTable path', type=str, default=None)
    parser.add_argument('--n_pcs', help='Number of PCs to be compute', type=int, default=20)
    parser.add_argument('--af_threshold', help='MAF cutoff for filtering', type=float, default=0.01)
    parser.add_argument('--ld_pruning', help='Perform LD pruning before PCA (recommended)', action='store_true')
    parser.add_argument('--r2', help='Squared correlation threshold for LD pruning', type=float, default=0.2)
    parser.add_argument('--bp_window_size', help='Window size in bps', type=int, default=500000)
    parser.add_argument('--export_loadings', help='Compute/export variant loadings as HailTable', action='store_true')
    parser.add_argument('--write_to_file', help='Write results to TSV file', action='store_true')

    args = parser.parse_args()

    main(args)
