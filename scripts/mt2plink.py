"""
Simple script to export PLINK files given a Hail MatrixTable.
Optionality, apply user-defined filters. Useful to run SKAT-O.

"""

import time
from config import HAIL_LOG_PATH
from utils.hail_functions import *
import hail as hl

# Initializing Hail on cluster mode
init_hail_on_cluster(tmp_dir='/mnt/nfs/mdatanode/hail-temp',
                     log_file=HAIL_LOG_PATH,
                     local_mode=False)

# Import MatrixTable
MT_INPUT_PATH = '/mnt/nfs/mdatanode/wes10k_resources/wes1k/mts/mt_1k_qced_filtered_09-07-2019.mt'
mt = hl.read_matrix_table(MT_INPUT_PATH)

# Filter MatrixTable before exporting (e.g. by AF)
rare_threshold = 0.01  # MAF 1%
mt = mt.filter_rows((mt.variant_qc.AF[1] <= rare_threshold) &
                    (mt['gnomad2_popmax_gnomad.AF'] <= rare_threshold) &
                    (mt['exac_af_adj'] <= rare_threshold) &
                    (mt.symbol == 'HEY2') &
                    (mt.csq_type != 'SYN'),
                    keep=True)

# Annotate sex info
mt = (mt
      .annotate_cols(is_female=hl.case()
                     .when(mt.f_stat <= 0.4, True)
                     .when(mt.f_stat >= 0.6, False)
                     .or_missing())
      )

# Export plink files
date = time.strftime("%d-%m-%Y")
PLINK_OUTPUT_PATH = f'/mnt/nfs/mdatanode/wes10k_resources/wes1k/plink_output/hey2_{date}'

hl.export_plink(dataset=mt,
                output=PLINK_OUTPUT_PATH,
                ind_id=mt.s,
                pheno=mt.isCase,
                is_female=mt.is_female)

# Export useful info (e.g. covariates, annotation)
delimiter = '|'
sample_expr_annotations = dict(
    cases_het=hl.delimit(hl.agg.filter(mt.GT.is_het() & mt.isCase, hl.agg.collect(mt.s)), delimiter=delimiter),
    cases_hom=hl.delimit(hl.agg.filter(mt.GT.is_hom_var() & mt.isCase, hl.agg.collect(mt.s)), delimiter=delimiter),
    controls_het=hl.delimit(hl.agg.filter(mt.GT.is_het() & ~mt.isCase, hl.agg.collect(mt.s)), delimiter=delimiter),
    controls_hom=hl.delimit(hl.agg.filter(mt.GT.is_hom_var() & ~mt.isCase, hl.agg.collect(mt.s)), delimiter=delimiter))

variant_table = (mt
                 .annotate_rows(**sample_expr_annotations)
                 .rows()
                 .select('key_variant',
                         'symbol',
                         'consequences',
                         'csq_type',
                         'cadd_phred',
                         'cases_het',
                         'cases_hom',
                         'controls_het',
                         'controls_hom')
                 )

# Export table
variant_table.export(output=f'{PLINK_OUTPUT_PATH}/variant_table.tsv')

# Stop Hail
hl.stop()
