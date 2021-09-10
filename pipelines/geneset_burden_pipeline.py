import argparse
import hail as hl
from config import HAIL_TMP_DIR, HAIL_LOG_PATH
from utils.hail_functions import init_hail_on_cluster, logistic_regression


def main(args):
    # Initializing Hail on cluster mode
    init_hail_on_cluster(tmp_dir=HAIL_TMP_DIR,
                         log_file=HAIL_LOG_PATH,
                         local_mode=True)

    # 1- Aggregate MatrixTable per gene/consequences creating gene/csq X sample matrix

    # Read MatrixTable
    mt = hl.read_matrix_table(args.mt_input_path)

    # Annotate csq group info per variants
    # Define consequences variant rules with hail expressions
    # TODO: check if field exist in dataset
    csq_group_rules = {}
    if args.ptv:
        csq_group_rules.update({'PTV': mt.csq_type == 'PTV'})
    if args.pav:
        csq_group_rules.update({'PAV': mt.csq_type == 'PAV'})
    if args.syn:
        csq_group_rules.update({'SYN': mt.csq_type == 'SYN'})
    if args.cadd:
        sq_group_rules.update({'CADD': (mt.csq_type == 'PAV') &
                                                          (mt.cadd_phred >= args.cadd_threshold)})
    if args.mpc:
        csq_group_rules.update({'MPC': (mt.csq_type == 'PAV') &
                                                         (mt.mpc >= args.mpc_threshold)})

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
          )

    # Group mt by gene/csq_group.
    mt_grouped = (mt
                  .group_rows_by(mt.csq_group, mt.symbol)
                  .partition_hint(100)
                  .aggregate(n_het=hl.agg.count_where(mt.GT.is_het()))
                  )

    # 2- Annotate gene set information

    # Import/parsing gene cluster table
    clusters = hl.import_table(args.gene_set_path,
                               no_header=True)

    # parsing gene set column
    clusters = (clusters
                .transmute(genes=hl.set(clusters['f1'].split(delim='[|]')))
                )

    clusters = (clusters
                .explode(clusters.genes)
                )

    clusters = (clusters
                .group_by('genes')
                .partition_hint(100)
                .aggregate(cluster_name=hl.agg.collect_as_set(clusters['f0']))
                .key_by('genes')
                )

    # annotate gene set info
    mt_grouped = (mt_grouped
                  .annotate_rows(cluster_name=clusters[mt_grouped.symbol].cluster_name)
                  )

    # 3- Aggregate per gene set and consequences

    # Group mt by gene set/csq_group.
    mt_grouped = (mt_grouped
                  .explode_rows(mt_grouped.cluster_name)
                  )
    mt_grouped = (mt_grouped
                  .group_rows_by(mt_grouped.cluster_name, mt_grouped.csq_group)
                  .partition_hint(100)
                  .aggregate(n_het=hl.agg.sum(mt_grouped.n_het))
                  )

    # force to eval all aggregation operation by writing mt to disk
    mt_grouped = mt_grouped.persist(storage_level='DISK_ONLY')

    if args.logistic_regression:
        # covariates list
        covs = list(args.covs_list)

        # Define x expression (entries/genotype)
        x_expr = 'n_het'

        extra_annotations = {'analysis': 'all_cases',
                             'covariates': covs}

        tb_stats = logistic_regression(mt=mt_grouped,
                                       x_expr=x_expr,
                                       response=args.phenotype_field,
                                       covs=covs,
                                       pass_through=[],
                                       extra_fields=extra_annotations)
        # export table
        tb_stats.export(args.output_path)

    if args.fet:
        None  # TODO: implement Fisher Exact-based burden gene set test

    hl.stop()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    # input/output paths arguments
    parser.add_argument('--mt_input_path', help='Path to MatrixTable', type=str, default=None)
    parser.add_argument('--gene_set_path', help='Path to two-columns TSV gene set file', type=str,
                        default=None)
    parser.add_argument('--output_path', help='Output path to TSV file', type=str,
                        default=None)

    # consequences groups
    parser.add_argument('--ptv', help='Test protein truncating variants', action='store_true')
    parser.add_argument('--pav', help='Test protein altering variant', action='store_true')
    parser.add_argument('--syn', help='Test synonymous variants', action='store_true')
    parser.add_argument('--cadd', help='Test PAVs filtered by CADD-score', action='store_true')
    parser.add_argument('--cadd_threshold', help='CADD-score lower threshold', type=float, default=20)
    parser.add_argument('--mpc', help='Test MPC filtered by MPC-score', action='store_true')
    parser.add_argument('--mpc_threshold', help='MPC-score lower threshold', type=float, default=2)

    # statistical test to run
    parser.add_argument('--logistic_regression', help='Run logistic regression burden test', action='store_true')
    parser.add_argument('--phenotype_field', help='Binary phenotype field name', type=str, default='isCase')
    parser.add_argument('--add_covariates', help='Run logistic regression test with covariates', action='store_true')
    parser.add_argument('--covs_list', help='List of covariates field names to run logistic regression test',
                        type=str, nargs="*", default=[])
    parser.add_argument('--fet', help='Run tow-sided Fisher Exact Test', action='store_true')

    # x-expression
    parser.add_argument('--only_het', help='Aggregate only heterozygous genotypes', action='store_true')

    args = parser.parse_args()

    main(args)
