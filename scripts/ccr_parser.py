
import time
import hail as hl
from utils.hail_functions import add_global_annotations
from utils.parser import get_files_names
from config import ROOT_DIR, HAIL_LOG_PATH
import os

"""

Simple script to parse BED-formatted Constrain Coding Regions (CCRs)
and export the it as Hail Table with global annotations.

"""

# Define files paths
CCR_FILES_PATH = os.path.join(ROOT_DIR, 'testdata/ccrs')
OUTPUT_PATH = os.path.join(ROOT_DIR, 'testdata/hailtables/ccrs_test_table.ht')

# Initializing Hail on local mode
hl.init(log=HAIL_LOG_PATH)

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
annotations = {'date': time.strftime("%d-%m-%Y %H:%M"),
               'table_info': 'constrained coding regions (ccr)',
               'DOI': 'https://doi.org/10.1038/s41588-018-0294-6',
               'source_bed_files': 'https://github.com/quinlan-lab/ccrhtml',
               'field_info': {'interval': 'locus interval',
                              'ccr_pct': 'ccr percentile'},
               'genome_reference': 'hg19'}

ccr_tb = add_global_annotations(ccr_tb, annotations=annotations, overwrite=True)

# Write table to disk
ccr_tb.write(OUTPUT_PATH, overwrite=True)

hl.stop()
