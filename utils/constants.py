
# Consequences (VEP API 94)
# Define major consequences groups

CSQ_PTV = ['frameshift_variant',
           'stop_gained',
           'splice_donor_variant',
           'splice_acceptor_variant']

CSQ_PAV = ['missense_variant',
           'inframe_deletion',
           'inframe_insertion',
           'start_lost',
           'stop_retained_variant',
           'stop_lost',
           'protein_altering_variant',
           'start_retained_variant']

CSQ_SYN = ['synonymous_variant']

CSQ_CODING = CSQ_PTV + CSQ_PAV + CSQ_SYN

CSQ_NON_CODING = ['non_coding_transcript_exon_variant',
                  'intron_variant',
                  'splice_region_variant',
                  'upstream_gene_variant',
                  'downstream_gene_variant',
                  '3_prime_UTR_variant',
                  '5_prime_UTR_variant',
                  'mature_miRNA_variant',
                  'intergenic_variant',
                  'non_coding_transcript_variant',
                  'coding_sequence_variant',
                  'incomplete_terminal_codon_variant']
