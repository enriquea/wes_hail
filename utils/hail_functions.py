"""

A set of helper function for Hail-based pipelines.

"""

import os
import sys
import time
from typing import *

import hail as hl
import hail.expr.aggregators as agg
import pyspark
from pyspark import SparkContext


# run MT rows logistic regression
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


# run two-tailed fisher exact test
def compute_fisher_exact(tb: hl.Table,
                         n_cases_col: str,
                         n_control_col: str,
                         total_cases: int,
                         total_control: int,
                         correct_total_counts: bool,
                         root_col_name: str,
                         extra_fields: dict) -> hl.Table:
    """
    Perform two-sided Fisher Exact test. Add extra annotations (if any)

    :param tb: Hail Table
    :param n_cases_col: number of cases
    :param n_control_col: number of control
    :param total_cases: total cases
    :param total_control: total controls
    :param correct_total_counts: should the total numbers (case/control) be corrected to avoid duplicated counting?
    :param root_col_name: field to be annotated with test results
    :param extra_fields: Extra filed (should be a dict) to be annotated
    :return: Hail Table with Fisher Exact test results.
    """
    # compute fisher exact
    if correct_total_counts:
        fet = hl.fisher_exact_test(c1=hl.int32(tb[n_cases_col]),
                                   c2=hl.int32(tb[n_control_col]),
                                   c3=hl.int32(total_cases) - hl.int32(tb[n_cases_col]),
                                   c4=hl.int32(total_control) - hl.int32(tb[n_control_col]))
    else:
        fet = hl.fisher_exact_test(c1=hl.int32(tb[n_cases_col]),
                                   c2=hl.int32(tb[n_control_col]),
                                   c3=hl.int32(total_cases),
                                   c4=hl.int32(total_control))

    tb = (tb
          .annotate(**{root_col_name: fet})
          .flatten()
          )

    if len(extra_fields) == 0:
        return tb
    else:
        return tb.annotate(**extra_fields)


# Add global annotations to HailTable/MatrixTable
def add_global_annotations(ds: Union[hl.MatrixTable, hl.Table],
                           annotations: dict,
                           overwrite: True) -> Union[hl.MatrixTable, hl.Table]:
    """
    Rewrite MatrixTable/Table global annotations

    :param ds: MatrixTable or HailTable
    :param annotations: Dictionary with key:value annotations
    :param overwrite: If True (default), old annotation will be removed
    :return: ds (input) with global annotations
    """
    if overwrite:
        ds = (ds
              .drop(*ds.index_globals())  # remove global annotation
              .annotate_globals(**annotations)
              )
    else:
        ds = (ds
              .annotate_globals(**annotations)
              )
    return ds


# Function to compute sample/variant QC
def compute_qc(mt: hl.MatrixTable,
               root_col_name='sample_qc',
               root_row_name='variant_qc') -> hl.MatrixTable:
    """
    Given a MatrixTable, compute samples/variants quality controls metrics

    :param mt: Input MatrixTable
    :param root_col_name: prefix sample qc field
    :param root_row_name: prefix variant qc field
    :return: MatrixTable with quality control computed
    """
    mt = hl.sample_qc(mt, name=root_col_name)
    mt = hl.variant_qc(mt, name=root_row_name)
    return mt


# Basic function to initialize Hail on cluster or local mode
def init_hail_on_cluster(app_name: str = 'Hail',
                         genome_ref: str = 'GRCh38',
                         tmp_dir: str = '/tmp',
                         log_file: str = 'logs/hail.log',
                         local_mode: bool = False) -> None:
    """
    Init Hail context on Spark cluster or in local mode

    :param genome_ref: object
    :param app_name: Application name
    :param tmp_dir: Directory for temporal file (network visible on cluster mode)
    :param log_file: Hail log file
    :param local_mode: Init Hail on local mode
    :return: None
    """

    if local_mode:
        # Init Hail with a basic SparkContext
        number_cores = 4
        memory_gb = 16
        conf = (pyspark.SparkConf()
                .setMaster('local[{}]'.format(number_cores))
                .set('spark.driver.memory', '{}g'.format(memory_gb))
                )
        sc = pyspark.SparkContext(conf=conf)
        hl.init(sc=sc, default_reference=genome_ref)
    else:
        # Create SparkContext with default parameters (from SPARK_PATH/conf/spark-defaults.conf)
        sc = SparkContext(appName=app_name)
        hl.init(sc=sc,
                default_reference=genome_ref,
                tmp_dir=tmp_dir,
                app_name=app_name,
                log=log_file)


# Annotate Coding Constraint Regions (CCRs)
def annotate_ccr(mt: hl.MatrixTable,
                 ht_ccr: hl.Table,
                 fields_to_annotate: List[str]) -> hl.MatrixTable:
    """
    Annotate a given MatrixTable with Coding Constrain Region (CCRs) information.

    :param mt: Input MatrixTable
    :param ht_ccr: CCR Hail Table
    :param fields_to_annotate: List of fields from CCR table to be annotated
    :return: Annotated MatrixTable
    """
    mt = mt.annotate_rows(**ht_ccr.select(*fields_to_annotate)[mt.locus])
    return mt


# import list of vcf files (paths)
def get_files_names(path: str, ext: str) -> List[str]:
    """
    Get full path for files (e.g. vcf) in a given directory

    :param path: directory global path
    :param ext: files extension (e.g. vcf)
    :return: List of files in directory
    """
    list_files: List[str] = []
    for (dir_name, _, files) in os.walk(path):
        for filename in files:
            if filename.endswith(ext):
                list_files.append(os.path.abspath(os.path.join(dir_name, filename)))
    return list_files


# Summarize Coding Constraint Region (CCRs) table
def summary_ccr(ht_ccr: hl.Table,
                file_output: str,
                ccr_pct_start: int = 0,
                ccr_pct_end: int = 100,
                ccr_pct_bins: int = 10,
                cumulative_histogram: bool = False,
                ccr_pct_cutoffs=None) -> None:
    """
    Summarize Coding Constrain Region information (as histogram) per gene.

    :param ht_ccr: CCR Hail table
    :param file_output: File output path
    :param ccr_pct_start: Start of histogram range.
    :param ccr_pct_end: End of histogram range
    :param ccr_pct_bins: Number of bins
    :param cumulative_histogram: Generate a cumulative histogram (rather than to use bins)
    :param ccr_pct_cutoffs: Cut-offs used to generate the cumulative histogram
    :return: None
    """

    if ccr_pct_cutoffs is None:
        ccr_pct_cutoffs = [90, 95, 99]

    if cumulative_histogram:
        # generate cumulative counts histogram
        summary_tb = (ht_ccr
                      .group_by('gene')
                      .aggregate(**{'ccr_above_' + str(ccr_pct_cutoffs[k]): agg.filter(ht_ccr.ccr_pct >=
                                                                                       ccr_pct_cutoffs[k], agg.count())
                                    for k in range(0, len(ccr_pct_cutoffs))})
                      )
    else:
        summary_tb = (ht_ccr
                      .group_by('gene')
                      .aggregate(ccr_bins=agg.hist(ht_ccr.ccr_pct, ccr_pct_start, ccr_pct_end, ccr_pct_bins))
                      )

        # get bin edges as list (expected n_bins + 1)
        bin_edges = summary_tb.aggregate(agg.take(summary_tb.ccr_bins.bin_edges, 1))[0]

        # unpack array structure and annotate as individual fields
        summary_tb = (summary_tb
                      .annotate(**{'ccr_bin_' + str(bin_edges[k]) + '_' + str(bin_edges[k + 1]):
                                       summary_tb.ccr_bins.bin_freq[k] for k in range(0, len(bin_edges) - 1)})
                      .flatten()
                      )

        # drop fields
        fields_to_drop = ['ccr_bins.bin_edges', 'ccr_bins.bin_freq']
        summary_tb = (summary_tb
                      .drop(*fields_to_drop)
                      )

    # Export summarized table
    (summary_tb
     .export(output=file_output)
     )


# Annotate sex based on inbreeding coefficient (F_stat) on chromosome X
def annotate_sex(ds: Union[hl.MatrixTable, hl.Table],
                 f_stat_field: str = 'f_stat',
                 sex_field: str = 'sex',
                 female_upper_threshold: float = 0.4,
                 male_lower_threshold: float = 0.6,
                 ) -> Union[hl.MatrixTable, hl.Table]:
    """
    Annotate sex (0-female/1-male) based on F_stat (inbreeding coefficient computed on chr X)

    :param ds: Input MatrixTable or HailTable
    :param f_stat_field: F stat field name
    :param sex_field: Sex field name to be annotated
    :param female_upper_threshold: F_stat female upper threshold
    :param male_lower_threshold: F_stat male lower threshold
    :return: Annotated ds
    """
    if isinstance(ds, hl.MatrixTable):
        ds = (ds
              .annotate_cols(**{sex_field: (hl.case()
                                            .when(ds[f_stat_field] <= female_upper_threshold, 0)
                                            .when(ds[f_stat_field] >= male_lower_threshold, 1)
                                            .or_missing())}
                             )
              )
    else:
        ds = (ds
              .annotate(**{sex_field: (hl.case()
                                       .when(ds[f_stat_field] <= female_upper_threshold, 0)
                                       .when(ds[f_stat_field] >= male_lower_threshold, 1)
                                       .or_missing())}
                        )
              )
    return ds


def hail_version() -> str:
    return hl.__version__


# Filter variants defined in target regions.
def filter_interval(mt: hl.MatrixTable,
                    tb_bed: hl.Table) -> hl.MatrixTable:
    """
    Keep variants defined in target regions (BED file)

    :param mt: Hail MatrixTable
    :param tb_bed: HailTable with interval field generated with hl.import_bed function
    :return: Row-filtered MatrixTable
    """
    n_total = mt.count_rows()

    mt = (mt
          .filter_rows(hl.is_defined((tb_bed[mt.locus])),
                       keep=True)
          )

    n_filtered = n_total - mt.count_rows()

    pct = round((n_filtered / n_total) * 100, 2)

    print(f"Filtered {n_filtered} ({pct}%) non-covered variants out of {n_total}")

    return mt


# Print current date
def current_date() -> str:
    return time.strftime("%d%m%Y")


# Annotate fields from array (all elements must have the same length)
def annotate_from_array(ht: hl.Table,
                        array_field: str,
                        field_names: list) -> hl.Table:
    """
    Expand an array structure and add new fields.

    :param ht: HailTable
    :param array_field: The array field to be expanded
    :param field_names: The pre-defined fields names (ordered). Number of fields must match with array length
    :return: Annotated HailTable
    """

    # number of fields to be annotated
    n_fields = field_names.__len__()

    # get array field length
    array_len = hl.len(ht[array_field]).take(1)[0]

    if array_len == n_fields:
        tb_expanded = ht.transmute(**{field_names[i]: ht[array_field][i] for i in range(n_fields)})
    else:
        print("Number of fields don't match with array length...")
        sys.exit()
    return tb_expanded


# Keep only bi-allelic variants
def filter_biallelic(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Remove multi-allelic variant. Keep bi-allelic only.

    :param mt: Hail MatrixTable
    :return: Filtered MatrixTable
    """
    n_total = mt.count_rows()
    mt = (mt
          .filter_rows(hl.len(mt.alleles) == 2,
                       keep=True)

          )
    n_filtered = n_total - mt.count_rows()

    pct = round((n_filtered / n_total) * 100, 2)

    print(f'Filtered {n_filtered} ({pct}%) multi-allelic variants out of {n_total}.')

    return mt


# annotate variant key
def annotate_variant_key(ds: Union[hl.MatrixTable, hl.Table]
                         ) -> Union[hl.MatrixTable, hl.Table]:
    # define key variant expression
    key_expr = hl.delimit([ds.locus.contig,
                           hl.str(ds.locus.position),
                           ds.alleles[0],
                           ds.alleles[1]], ':')

    if isinstance(ds, hl.MatrixTable):
        ds = ds.annotate_rows(variant_key=key_expr)

    if isinstance(ds, hl.Table):
        ds = ds.annotate(variant_key=key_expr)

    return ds
