import argparse
import hail as hl


# function to run MT rows logistic regression
def logistic_regression(mt: hl.MatrixTable,
                        x_expr: str,
                        response: str,
                        covs: list,
                        pass_through: list,
                        extra_fields: dict,
                        add_odd_stats: bool = True) -> hl.Table:
    """
    Perform a logistic-regression test (use by default Wald test).

    :param mt: Hail MatrixTable
    :param x_expr: the genotype field name (numeric expression)
    :param response: binary response
    :param covs: list of covariates to be included in the test
    :param pass_through: list of extra fields to keep in the output
    :param extra_fields: extra field to annotated (expected a dict)
    :param add_odd_stats: compute odds from logistic regression stats
    :return: Hail Table with logistic regression test results
    """
    # parsing covariates list
    if len(covs) >= 1:
        covs = [1] + [mt[field] for field in covs]
    else:
        covs = [1]

    tb_stats = hl.logistic_regression_rows(y=mt[response],
                                           x=mt[x_expr],
                                           covariates=covs,
                                           pass_through=pass_through,
                                           test='wald')

    if add_odd_stats:
        # Compute Odds ratio and 95% confidence interval from logistics regression stats
        tb_stats = tb_stats.annotate(odds_ratio=hl.exp(tb_stats.beta),
                                     lower_ci_95=hl.exp(tb_stats.beta - 1.96 * tb_stats.standard_error),
                                     upper_ci_95=hl.exp(tb_stats.beta + 1.96 * tb_stats.standard_error))

    # add column with additional information
    if len(extra_fields) == 0:
        return tb_stats
    else:
        return tb_stats.annotate(**extra_fields)


def main(args):

    # Initializing Hail on cluster mode
    hl.init()

    # 1- Aggregate MatrixTable per gene/consequences creating gene/csq X sample matrix
    # Read MatrixTable
    mt = hl.read_matrix_table(args.mt_input_path)

    # Annotate csq group info per variants
    # Define consequences variant rules with hail expressions
    # TODO: check if fields exist in dataset
    csq_group_rules = {}
    if args.ptv:
        csq_group_rules = csq_group_rules.update({'PTV': mt.csq_type == 'PTV'})
    if args.pav:
        csq_group_rules = csq_group_rules.update({'PAV': mt.csq_type == 'PAV'})
    if args.syn:
        csq_group_rules = csq_group_rules.update({'SYN': mt.csq_type == 'SYN'})
    if args.cadd:
        csq_group_rules = csq_group_rules.update({'CADD': (mt.csq_type == 'PAV') &
                                                          (mt.cadd_phred >= args.cadd_threshold)})

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

    # force to eval all aggregation operation by writing mt to disk
    # mt_grouped = mt_grouped.persist(storage_level='DISK_ONLY')

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
        None  # TODO: implement gene-based Fisher Exact burden test

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
