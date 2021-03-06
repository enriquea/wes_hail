"""
Compute relatedness estimates between individuals using a variant of the PC-Relate method.

Example usage:

python sample_relatedness --mt_path 'path/to/mt' \
                          --ht_output_path 'path/to/hts' \
                          --common_only \
                          --maf_threshold 0.05 \
                          --compute_pcs \
                          --ld_pruning \
                          --write_to_file
"""

import argparse
import time
from config import *
from utils.hail_functions import *


def main(args):
    # Initializing Hail on cluster mode
    init_hail_on_cluster(app_name='sample-relatedness',
                         tmp_dir=HAIL_TMP_DIR,
                         log_file=HAIL_LOG_PATH,
                         local_mode=False)

    # Read MatrixTable
    mt = hl.read_matrix_table(args.mt_path)

    # force write to disk with less number of partition
    mt = mt.repartition(100, shuffle=False)
    mt.write('/mnt/nfs/mdatanode/hail-temp/temp.mt', overwrite=True)
    mt = hl.read_matrix_table('/mnt/nfs/mdatanode/hail-temp/temp.mt')

    if args.common_only:
        # Filter by AF (i.e. keep common variants)
        mt = (mt
              .filter_rows(mt.variant_qc.AF[1] >= args.maf_threshold, keep=True)
              )

    if args.sample_to_keep is not None:
        sample_table = hl.import_table(paths=args.sample_to_keep,
                                       no_header=True)
        sample_set = hl.literal(sample_table.aggregate(agg.collect_as_set(sample_table.f0)))
        mt = mt.filter_cols(sample_set.contains(mt.s), keep=True)

    if args.ld_pruning:
        # LD pruning
        # Avoid filtering / missingness entries (genotypes) before run LP pruning
        # Zulip Hail support issue -> "BlockMatrix trouble when running pc_relate"
        mt = mt.unfilter_entries()

        # Prune variants in linkage disequilibrium.
        # Return a table with nearly uncorrelated variants
        pruned_variant_table = hl.ld_prune(mt.GT,
                                           r2=args.r2,
                                           bp_window_size=args.bp_window_size,
                                           memory_per_core=512)

        # Keep LD-pruned variants
        mt = (mt
              .filter_rows(hl.is_defined(pruned_variant_table[mt.locus, mt.alleles]), keep=True)
              )

    relatedness = hl.pc_relate(mt.GT,
                               args.individual_specific_maf,
                               k=args.n_pcs,
                               statistics='all')

    # TODO: retrieve maximal independent sample set

    # Write/export table to file
    tb = (relatedness
          .flatten()
          .key_by('i.s', 'j.s')
          )
    date = time.strftime("%d%m%Y")
    tb.write(output=f'{args.ht_output_path}_{date}.ht')

    # Write PCs table to file (if specified)
    if args.write_to_file:
        # Export table to file
        tb.export(output=f'{args.ht_output_path}_kinship_stats_{date}.tsv')

    hl.stop()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--mt_path', help='Path to MatrixTable with computed QC metrics', type=str, default=None)
    parser.add_argument('--ht_output_path', help='Output HailTable path with Kinship stats', type=str, default=None)
    parser.add_argument('--sample_to_keep', help='Text file (one-column, no header) listing the samples to keep',
                        type=str, default=None)
    parser.add_argument('--n_pcs', help='Number of PCs to be computed', type=int, default=10)
    parser.add_argument('--individual_specific_maf', help='MAF cutoff for filtering', type=float, default=0.01)
    parser.add_argument('--common_only', help='Run PC-relate method on common variant set (default MAF 1%)',
                        action='store_true')
    parser.add_argument('--maf_threshold', help='Variants with maf below that threshold are removed',
                        type=float, default=0.01)
    parser.add_argument('--ld_pruning', help='Perform LD pruning before PCA (recommended)', action='store_true')
    parser.add_argument('--r2', help='Squared correlation threshold for LD pruning', type=float, default=0.2)
    parser.add_argument('--bp_window_size', help='Window size in bps', type=int, default=500000)
    parser.add_argument('--write_to_file', help='Write results to TSV file', action='store_true')

    args = parser.parse_args()

    main(args)
