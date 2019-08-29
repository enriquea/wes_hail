import argparse
from config import HAIL_LOG_PATH
from utils.hail_functions import *


def main(args):

    # Initializing Hail on cluster mode
    init_hail_on_cluster(tmp_dir='/mnt/nfs/mdatanode/hail-temp',
                         log_file=HAIL_LOG_PATH,
                         local_mode=False)

    # Read MatrixTable
    mt = hl.read_matrix_table(args.mt_input_path)

    # Import/parsing gene cluster table
    clusters = hl.import_table(args.gene_set_path)

    # parsing gene set column
    clusters = (clusters
                .transmute(gene_set=hl.set(clusters.geneset.split(delim='[|]')))
                .explode(clusters.gene_set)
                .group_by('gene_set')
                .partition_hint(100)
                .aggregate(cluster_name=agg.collect_as_set(clusters.cluster_name))
                .key_by('gene_set')
                )

    # annotate gene set info
    mt = (mt
          .annotate_rows(gene_set=clusters[mt.symbol].cluster_name)
          )

    # Annotate csq group info per variants
    # Define consequences variant rules with hail expressions
    csq_group_rules = {'PTV': mt.csq_type == 'PTV',
                       'PAV': mt.csq_type == 'PAV',
                       'SYN': mt.csq_type == 'SYN',
                       'CADD20': (mt.csq_type == 'PAV') & (mt.cadd_phred >= args.cadd_threshold),
                       'MPC2': (mt.csq_type == 'PAV') & (mt.mpc >= args.mpc_threshold)
                       }

    # Annotate groups per variants
    mt = (mt
          .annotate_rows(csq_group=csq_group_rules)
          )

    # Transmute csq_group and convert to set (easier to explode and grouping later)
    mt = (mt
          .transmute_rows(csq_group=hl.set(hl.filter(lambda x:
                                                     mt.csq_group.get(x),
                                                     mt.csq_group.keys())))
          )

    # Explode nested csq_group before grouping
    mt = (mt
          .explode_rows(mt.csq_group)
          .explode_rows(mt.cluster_name)
          )

    # Group mt by cluster id
    mt_grouped = (mt
                  .group_rows_by(mt.cluster_name, mt.csq_group)
                  .partition_hint(100)
                  .aggregate(n_het=agg.count_where(mt.GT.is_het()))
                  )

    if args.logistic_regression:
        # covariates list
        covs = list(args.covs_list)

        # Define x expression (entries/genotype)
        x_expr = 'n_het'

        # running syndromic
        extra_annotations = {'analysis': 'all_cases',
                             'covariates': covs}

        tb_stats = logistic_regression(mt=mt_grouped,
                                       x_expr=x_expr,
                                       response='isCase',
                                       covs=covs,
                                       pass_through=[],
                                       extra_fields=extra_annotations)
        # export table
        tb_stats.export(args.output_path)

    if args.fet:
        None  # TODO: implement Fisher Exact-based burden gene set test


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    # input/output paths arguments
    parser.add_argument('--mt_input_path', help='Path to MatrixTable', type=str, default=None)
    parser.add_argument('--gene_set_path', help='Path to two-columns TSV gene set file', type=str,
                        default=None)
    parser.add_argument('--output_path', help='Output path to TSV file', type=str,
                        default=None)

    # consequences groups
    # parser.add_argument('--ptv', help='Test protein truncating variants', action='store_true')
    # parser.add_argument('--pav', help='Test protein altering variant', action='store_true')
    # parser.add_argument('--syn', help='Test synonymous variants', action='store_true')
    parser.add_argument('--cadd', help='Test PAVs filtered by CADD-score', action='store_true')
    parser.add_argument('--cadd_threshold', help='CADD-score lower threshold', type=float, default=20)
    parser.add_argument('--mpc', help='Test MPC filtered by MPC-score', action='store_true')
    parser.add_argument('--mpc_threshold', help='MPC-score lower threshold', type=float, default=20)

    # test to run
    parser.add_argument('--logistic_regression', help='Run logistic regression burden test', action='store_true')
    parser.add_argument('--phenotype_field', help='Binary phenotype field name', type=str, default='isCase')
    parser.add_argument('--add_covariates', help='Run logistic regression test with covariates', action=True)
    parser.add_argument('--covs_list', help='List of covariates field names to run logistic regression test', type=list,
                        default=None)
    parser.add_argument('--fet', help='Run tow-sided Fisher Exact Test', action='store_true')

    # x-expression
    parser.add_argument('--only_het', help='Aggregate only heterozygous genotypes', action=True)

    args = parser.parse_args()

    main(args)
