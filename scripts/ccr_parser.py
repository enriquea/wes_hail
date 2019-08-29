"""

Simple script to parse BED-formatted Constrain Coding Regions (CCRs)
and export the it as Hail Table with global annotations.

"""

import time

from config import HAIL_LOG_PATH
from utils.hail_functions import *

# Define files paths
CCR_FILES_PATH = '/mnt/nfs/mdatanode/wes10k_resources/ccrs'
OUTPUT_PATH = '/mnt/nfs/mdatanode/wes10k_resources/ccrs/ccr_table_072019.ht'

# Initializing Hail on cluster mode
init_hail_on_cluster(tmp_dir='/mnt/nfs/mdatanode/hail-temp',
                     log_file=HAIL_LOG_PATH,
                     local_mode=False)

# Read CCRs BED files
list_ccr_files = get_files_names(path=CCR_FILES_PATH,
                                 ext='bed.gz')

# Rename chromosome field name
ccr_tb = (hl.import_table(list_ccr_files, force_bgz=True)
          .rename({'#chrom': 'chr'})
          )

# Coerce fields to expected types
ccr_tb = (ccr_tb
          .transmute(chr=hl.str(ccr_tb.chr),
                     start=hl.int(ccr_tb.start),
                     end=hl.int(ccr_tb.end),
                     ccr_pct=hl.float(ccr_tb.ccr_pct))
          )

# Making locus intervals from contig and positions
ccr_tb = (ccr_tb
          .annotate(interval=hl.parse_locus_interval(hl.str(ccr_tb.chr +
                                                            ':' + hl.str(ccr_tb.start) +
                                                            '-' + hl.str(ccr_tb.end))))
          .key_by('interval')
          )

# Annotate global annotations
annotations = {'date': time.strftime("%d-%m-%Y"),
               'table_info': 'constrained coding regions (ccr)',
               'DOI': 'https://doi.org/10.1038/s41588-018-0294-6',
               'source_bed_files': 'https://github.com/quinlan-lab/ccrhtml',
               'field_info': {'interval': 'locus interval',
                              'ccr_pct': 'ccr percentile'},
               'genome_reference': 'hg19'}

ccr_tb = add_global_annotations(ccr_tb, annotations=annotations, overwrite=True)

# Print table fields
ccr_tb.describe()

# Write table to disk
ccr_tb.write(OUTPUT_PATH, overwrite=True)

# Stop cluster
hl.stop()
