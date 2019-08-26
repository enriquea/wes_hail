"""

A set of helper function for Hail-based pipelines.

"""

import os
from typing import *
from typing import List

import hail as hl
import hail.expr.aggregators as agg
from pyspark import SparkContext


# run MT rows logistic regression
def logistic_regression(mt: hl.MatrixTable,
                        x_expr: str,
                        response: str,
                        covs: list,
                        pass_through: list,
                        extra_fields: dict,
                        add_odd_stats: bool = True) -> hl.Table:
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
                         root_col_name: str,
                         corrected_total_count: True,
                         extra_fields: dict) -> hl.Table:
    # compute fisher exact
    if corrected_total_count:
        fet = hl.fisher_exact_test(c1=hl.int32(tb[n_cases_col]),
                                   c2=hl.int32(tb[n_control_col]),
                                   c3=hl.int32(total_cases - tb[n_cases_col]),
                                   c4=hl.int32(total_control - tb[n_control_col]))
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
                         tmp_dir: str = '/tmp',
                         log_file: str = '/logs/hail.log',
                         local_mode: bool = False):
    if local_mode:
        hl.init()
    else:
        # Create SparkContext with default parameters (from SPARK_PATH/conf/spark-defaults.conf)
        sc = SparkContext(appName=app_name)
        hl.init(sc=sc,
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
