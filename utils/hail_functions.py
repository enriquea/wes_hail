import hail as hl
import hail.expr.aggregators as agg
from typing import Union


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
