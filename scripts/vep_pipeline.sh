#!/usr/bin/env bash

set -e
set -u
set -o pipefail

# Date: 27.02.2020
# This script holds the routines to run the VEP annotation tool.
# It adds some custom annotations and uses different plugins (e.g. CADD).

# Option used description here (--everything):
# --sift b, --polyphen b, --ccds, --uniprot, --hgvs, --symbol, --numbers,
# --domains, --regulatory, --canonical, --protein, --biotype, --uniprot, --tsl,
# --appris, --gene_phenotype --af, --af_1kg, --af_esp, --af_gnomad, --max_af,
# --pubmed, --variant_class, --mane

# change to dir with VCF files
cd ~/projects/github/wes_hail/testdata/

# set path for cache, plugins, custom annotation sources, etc...
vep_dir=~/.vep # path to ensembl cache data
custom_dir=/media/biolinux/BData/vep_plugins_cache # path to custom annotation files

for inputfile in /media/biolinux/BData/wes_ukbb_chd_50k/*_vcf.bgz
 do
  outfile=$(basename "$inputfile")
  outfile="chr21_veped_${outfile}"

  echo Processing ${inputfile} file...

  perl ~/tools/ensembl99-vep/vep -i "$inputfile" -o ${outfile} \
       --assembly 'GRCh38' \
       --dir ${vep_dir} \
       --fork 4 \
       --buffer_size 5000 \
       --cache \
       --cache_version 99 \
       --offline \
       --everything \
       --plugin CADD,${custom_dir}/cadd/hg38/whole_genome_SNVs.tsv.gz,${custom_dir}/cadd/hg38/InDels.tsv.gz \
       --custom ${custom_dir}/clinvar/hg38/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG \
       --custom ${custom_dir}/gnomad/hg38/gnomad.genomes.r3.0.sites.chr21.vcf.bgz,gnomADg,vcf,exact,0,AF,AC,AF_afr,AC_afr,AF_ami,AC_ami,AF_amr,AC_amr,AF_asj,AC_asj,AF_eas,AC_eas,AF_fin,AC_fin,AF_nfe,AC_nfe,AF_sas,AC_sas,AF_oth,AC_oth \
       --vcf \
       --compress_output 'bgzip' \
       --stats_text \
       --force_overwrite
 done