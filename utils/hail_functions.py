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
                        extra_fields: dict) -> hl.Table:
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
