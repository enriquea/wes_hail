# Script to extract and filtering the 1K CHD-cohort for replication

# Start a new Hail context

import time

from config import HAIL_LOG_PATH
from utils.hail_functions import *

# Initializing Hail on cluster mode
init_hail_on_cluster(tmp_dir='/mnt/nfs/mdatanode/hail-temp',
                     log_file=HAIL_LOG_PATH,
                     local_mode=False)

# Read full annotated and pre-filtered MatrixTable
PATH_INPUT_MT = '/mnt/nfs/mdatanode/wes10k_resources/hail/mts/wes10k_pre_filtered_qced_14112018.mt'
mt = hl.read_matrix_table(PATH_INPUT_MT)

# Getting 1K subset (cases-control) for replication
# samples to keep
# Import TSV file as Hail Table
PATH_SAMPLE_FILE = '/mnt/nfs/mdatanode/wes10k_resources/wes1k/annotations/samples_ids_1k_study.tsv'
sample_tb = (hl.import_table(PATH_SAMPLE_FILE,
                             impute=True,
                             no_header=False)
             .key_by('ms_id'))

# annotate 1k set info
mt = mt.annotate_cols(**sample_tb[mt.s])

# collect samples as set
# sample_ids = 'MS_ID' # sample ids column name
# sample_set = hl.literal([row[sample_ids] for row in sample_tb.select(sample_ids).collect()])

# Keep only 1K subset
mt_1k = mt.filter_cols(hl.literal({'1K', 'Control'}).contains(mt['case_control']), keep=True)

# Re-compute QC based on this subset (e.g. Internal allelic frequency, Depth mean, Genotype quality mean, ect.)
mt_1k = compute_qc(mt_1k)

# Filtering by QC

# Filtering Samples
mt_1k = mt_1k.filter_cols((mt_1k.sample_qc.dp_stats.mean >= 15) &  # depth
                          (mt_1k.sample_qc.call_rate >= 0.95) &  # call rate
                          (mt_1k.sample_qc.gq_stats.mean >= 25) &  # quality
                          (mt_1k.max_pi_hat < 0.8) &  # relatedness
                          (mt_1k.probability_EUR > 0.65),  # ancestry
                          keep=True)

# Filtering by variants
mt_1k = (mt_1k
         .filter_rows((mt_1k.variant_qc.dp_stats.mean >= 15) &
                      (mt_1k.variant_qc.call_rate >= 0.95) &
                      (mt_1k.variant_qc.gq_stats.mean >= 25),
                      keep=True))

# Filtering by Genotype
ab = mt_1k.AD[1] / hl.sum(mt_1k.AD)

filter_condition_ab = ((mt_1k.GT.is_hom_ref() & (ab <= 0.1)) |
                       (mt_1k.GT.is_het() & (ab >= 0.25) & (ab <= 0.75)) |
                       (mt_1k.GT.is_hom_var() & (ab >= 0.9)))

mt_1k = mt_1k.filter_entries(filter_condition_ab, keep=True)

# Filtering by Rarity criteria
# rare_threshold = 0.01  # MAF 1%
# mt_rare = mt_1k.filter_rows((mt_1k.variant_qc.AF[1] <= rare_threshold), keep=True)

# Remove old global annotations and add new ones with the QC/filtering steps

ds = mt_1k

# cases/control counts
total_cases = ds.aggregate_cols(agg.count_where(ds.isCase))
total_control = ds.aggregate_cols(agg.count_where(~ds.isCase))
total_syndromic = ds.aggregate_cols(agg.count_where(ds.SyndromicCHD))
total_nonsyndromic = ds.aggregate_cols(agg.count_where(ds.NonSyndromicCHD))

# define global annotation
date = time.strftime("%d-%m-%Y")
annotations = {'date': date,
               'vep_api': 'v90',
               'cohort': 'wes1k-replication',
               'sample_pre_filter': 'removed duplicated 0.8|only Europen 0.65',
               'variant_pre_filter': 'covered variant (v3/v4/v5)|biallelic|only exonic',
               'QC_computing': 'QC computed after pre-filtering',
               'sample_qc_filter': 'call_rate >= 0.95|dp_mean >= 15|gq_mean >= 25',
               'variant_qc_filter': 'call_rate >= 0.95|dp_mean >= 15|gq_mean >= 25',
               'genotype_qc_filter': {'hets': '0.25 <= ab <= 0.75',
                                      'homs_ref': 'ab <= 0.10',
                                      'homs_var': 'ab >= 0.90'},
               'AF_filtering': 'None',
               'n_cases': total_cases,
               'n_syndromic': total_syndromic,
               'n_nonsyndromic': total_nonsyndromic,
               'n_control': total_control,
               'group_description': {'PTV': 'protein truncating variant',
                                     'PAV': 'protein altering variant',
                                     'SYN': 'silent'}}
# rewrite global annotations
ds = add_global_annotations(ds, annotations=annotations, overwrite=True)

# export annotated and filtered MatrixTable for post-processing
PATH_OUTPUT_MT = f'/mnt/nfs/mdatanode/wes10k_resources/wes1k/mts/mt_1k_qced_filtered_{date}.mt'
(ds
 .write(PATH_OUTPUT_MT,
        overwrite=False))

# Stop Hail
hl.stop()
