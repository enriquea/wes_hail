"""

Given a QCed MatrixTable, apply a set of user-defined filters (samples/variants/genotypes)

"""

import argparse
import functools
import operator
import time

from utils.hail_functions import *


def main(args):

    # Initializing Hail on cluster mode
    init_hail_on_cluster(local_mode=True)

    # Read MatrixTable
    mt = hl.read_matrix_table(args.mt_input_path)

    # applying sample filters
    filter_sample_expr = [mt.sample_qc.dp_stats.mean >= args.min_sample_dp_mean,
                          mt.sample_qc.call_rate >= args.min_sample_callrate,
                          mt.sample_qc.gq_stats.mean >= args.min_sample_gq_mean]

    mt = (mt
          .filter_cols(functools.reduce(operator.iand, filter_sample_expr), keep=True)
          )

    # applying variant filters
    filter_variant_expr = []

    if args.only_snp:
        filter_variant_expr.append(hl.is_snp(mt.alleles[0], mt.alleles[1]))
    if args.only_biallelic:
        filter_variant_expr.append(hl.len(mt.alleles) == 2)

    filter_variant_expr.append(mt.variant_qc.dp_stats.mean >= args.min_variant_dp_mean)
    filter_variant_expr.append(mt.variant_qc.call_rate >= args.min_variant_callrate)
    filter_variant_expr.append(mt.variant_qc.gq_stats.mean >= args.min_variant_gq_mean)

    mt = (mt
          .filter_rows(functools.reduce(operator.iand, filter_variant_expr), keep=True)
          )

    # genotype filtering

    if args.keep_allele_balanced_genotypes:
        ab = mt.AD[1] / hl.sum(mt.AD)  # expression to compute allelic balance

        filter_condition_ab = [mt.GT.is_hom_ref() & (ab <= 0.1),
                               mt.GT.is_het() & (ab >= 0.25) & (ab <= 0.75),
                               mt.GT.is_hom_var() & (ab >= 0.9)]

        mt = (mt
              .filter_entries(functools.reduce(operator.ior, filter_condition_ab), keep=True)
              )

    # annotated applied filter before exporting filtered matrix

    # define global annotation
    date = time.strftime("%d-%m-%Y")
    annotations = {'date': date,
                   'vep_api': 'v90',
                   'sample_qc_filter': {'min_call_rate': args.min_sample_callrate,
                                        'min_gq_mean': args.min_sample_callrate,
                                        'min_dp_mean': args.min_sample_dp_mean},
                   'variant_qc_filter': {'min_call_rate': args.min_variant_callrate,
                                         'min_gq_mean': args.min_variant_callrate,
                                         'min_dp_mean': args.min_variant_dp_mean},
                   'only_snp': True if args.only_snp else False,
                   'only_biallelic': True if args.only_biallelic else False,
                   'genotype_qc_filter': {'filtered_genotypes_by_ab': args.keep_allele_balanced_genotypes},
                   'AF_filtering': 'None'}

    # rewrite global annotations
    mt = add_global_annotations(mt, annotations=annotations, overwrite=True)

    (mt
     .write(args.mt_output_path,
            overwrite=True)
     )

    # Stop Hail
    hl.stop()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    # input/output paths arguments
    parser.add_argument('--mt_input_path', help='Path to MatrixTable with computed QC metrics', type=str, default=None)
    parser.add_argument('--mt_output_path', help='Output path to annotated and filtered MatrixTable', type=str,
                        default=None)

    # sample filters thresholds
    parser.add_argument('--min_sample_dp_mean', help='Minimal sample depth mean', type=float, default=15)
    parser.add_argument('--min_sample_callrate', help='Minimal sample call rate', type=float, default=0.95)
    parser.add_argument('--min_sample_gq_mean', help='Minimal sample genotype quality mean', type=float, default=25)

    # variant filters thresholds
    parser.add_argument('--only_biallelic', help='Keep only biallelic variants', action='store_true')
    parser.add_argument('--only_snp', help='Keep only SNP variants (e.g. remove indels)', action='store_true')
    parser.add_argument('--min_variant_dp_mean', help='Minimal variant depth mean', type=float, default=15)
    parser.add_argument('--min_variant_callrate', help='Minimal variant call rate', type=float, default=0.95)
    parser.add_argument('--min_variant_gq_mean', help='Minimal variant genotype quality mean', type=float, default=25)

    # genotype filtering by allelic balance
    parser.add_argument('--keep_allele_balanced_genotypes', help="Filter out non-balanced genotypes",
                        action='store_true')

    args = parser.parse_args()

    main(args)
