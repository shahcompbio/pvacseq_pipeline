import os
import pandas as pd

cohort_file = config['cohort_file']
cohort = pd.read_table(cohort_file) # patient, sample, rna_sample w/o NA

if config['test_mode']:
    # Test
    cohort = pd.DataFrame([['test', 'testDNA', 'testRNA']], # test data
                   columns=['patient', 'sample', 'rna_sample'])
   
pvacseq_flt_tsvs = [f'main_run/{patient}/{sample}/outputs/pvacseq/MHC_Class_I/{sample}.expression_filtered.tsv' 
            for rix, (patient, sample, rna_sample) in cohort.iterrows()]
rule all:
    input:
        pvacseq_flt_tsvs

include: "rules/pvacseq.smk"
