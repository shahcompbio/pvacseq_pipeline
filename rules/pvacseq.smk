import os
import pandas as pd

def _get_rna_id(sample, cohort=cohort):
    samples = cohort[cohort['sample']==sample]['rna_sample']
    assert len(samples) == 1, f'samples = {samples}'
    return samples.values[0]

def _get_kallisto(wildcards):
    rna_sample_id = _get_rna_id(wildcards.sample)
    paths_path = "/juno/work/shah/users/chois7/apollo/APOLLO.KALLISTO.tsv"
    assert os.path.exists(paths_path)
    paths = pd.read_table(paths_path)
    paths = paths[paths['isabl_sample_id']==rna_sample_id]
    paths = paths[ (paths["result_type"] == "abundance_tsv")]['result_filepath'].values
    assert len(paths) == 1, f'abundance_tsv paths length not 1 for {rna_sample_id}: {paths}'
    path = paths[0]
    return path

def _get_bam(wildcards):
    paths_path = '/juno/work/shah/users/chois7/apollo/APOLLO.WGS-ALIGNMENT.tsv'
    assert os.path.exists(paths_path)
    paths = pd.read_table(paths_path)
    paths = paths[paths['isabl_sample_id']==wildcards.sample]
    paths = paths[ (paths["result_type"] == "bam")]['result_filepath'].values
    assert len(paths) == 1, f'bam paths length not 1: {paths}'
    path = paths[0]
    return path

def _get_bai(wildcards):
    paths_path = '/juno/work/shah/users/chois7/apollo/APOLLO.WGS-ALIGNMENT.tsv'
    assert os.path.exists(paths_path)
    paths = pd.read_table(paths_path)
    paths = paths[paths['isabl_sample_id']==wildcards.sample]
    paths = paths[ (paths["result_type"] == "bai")]['result_filepath'].values
    assert len(paths) == 1, f'bai paths length not 1: {paths}'
    path = paths[0]
    return path
    
def _get_maf(wildcards):
    paths_path = '/juno/work/shah/users/chois7/apollo/APOLLO.WGS-SOMATICCALLING.tsv'
    assert os.path.exists(paths_path)
    paths = pd.read_table(paths_path)
    paths = paths[paths['isabl_sample_id']==wildcards.sample]
    paths = paths[ (paths["result_type"] == "consensus_somatic_maf")]['result_filepath'].values
    assert len(paths) == 1, f'consensus_somatic_maf paths length not 1: {paths}'
    path = paths[0]
    return path

 # extract hla reads with samtools view and sort
rule get_hla_reads:
    input:
        bam=_get_bam, # tumor bam 
        bam_index=_get_bai,
    output: # main_run/A002/ADT002/outputs/get_hla_reads/hla.bam
        hla_bam='main_run/{patient}/{sample}/outputs/get_hla_reads/hla.bam',
        hla_bai='main_run/{patient}/{sample}/outputs/get_hla_reads/hla.bam.bai',
    params:
        intervals=config['gatk']['interval'],
        gatk_hla_bai='main_run/{patient}/{sample}/outputs/get_hla_reads/hla.bai',
    singularity:
        'docker://broadinstitute/gatk',
    shell:
        #'{gatk} PrintReads '.format(gatk=config['gatk']['run_file']) 
        'gatk PrintReads ' 
        '-I {input.bam} ' 
        '-L {params.interval} ' # chrO:0-9 format
        '-O {output.hla_bam} && ' 
        'mv {params.gatk_hla_bai} {output.hla_bai}' # fix bai name

# optitype samtools conversion step
rule optitype_convert:
    input:
        hla_bam='main_run/{patient}/{sample}/outputs/get_hla_reads/hla.bam',
        hla_bai='main_run/{patient}/{sample}/outputs/get_hla_reads/hla.bam.bai',
    output: # main_run/A002/ADT002/outputs/optitype_convert_all
        chr6=temp('main_run/{patient}/{sample}/outputs/optitype_convert_all/sample_chr6.fastq'),
        unmapped=temp('main_run/{patient}/{sample}/outputs/optitype_convert_all/sample_unmapped.fastq'),
        chr6_unmapped=temp('main_run/{patient}/{sample}/outputs/optitype_convert_all/sample_chr6_unmapped.fastq'),
        fastq_R1='main_run/{patient}/{sample}/outputs/optitype_convert_all/sample_R1.fastq',
        fastq_R2='main_run/{patient}/{sample}/outputs/optitype_convert_all/sample_R2.fastq',
    singularity:
        "docker://biocontainers/samtools:v1.9-4-deb_cv1",
    shell:
        # Subset reads mapped to chr 6
        'samtools view -b -f3 {input.hla_bam} \"6\" > {output.chr6}; ' 
        # Subset unmapped reads
        'samtools view -b -f13 {input.hla_bam} > {output.unmapped}; ' 
        # Merge reads mapped to chr 6 and unmapped reads
        '(samtools merge -f {output.chr6_unmapped} {output.chr6} {output.unmapped}); ' 
        # Splitting BAM file into paired-ends and create FASTQs
        'samtools view -bf 0x40 {output.chr6_unmapped} | samtools bam2fq - > {output.fastq_R1}; ' 
        'samtools view -bf 0x80 {output.chr6_unmapped} | samtools bam2fq - > {output.fastq_R2}'

# optitype razers
rule optitype_extract:
    input:
        fastq_R1='main_run/{patient}/{sample}/outputs/optitype_convert_all/sample_R1.fastq',
        fastq_R2='main_run/{patient}/{sample}/outputs/optitype_convert_all/sample_R2.fastq',
    output: # main_run/A002/ADT002/outputs/optitype_razers/sample_R1.bam
        hla_bam_R1='main_run/{patient}/{sample}/outputs/optitype_razers/sample_R1.bam',
        hla_bam_R2='main_run/{patient}/{sample}/outputs/optitype_razers/sample_R2.bam',
    params:
        name='optitype-razers',
        perc_indent=config['hla']['optitype']['perc_indent'],
        max_hits=config['hla']['optitype']['max_hits'],
        dist_range=config['hla']['optitype']['dist_range'],
        hla_ref=config['hla']['optitype']['hla_ref'],
    singularity:
        "docker://fred2/optitype",
    shell:
        # create the command line using the parameters contained in config['razers']
        # add the output, reference and input paths that Snakemake will automatically replace
        '/usr/local/bin/razers3 -i {params.perc_indent} -m {params.max_hits} -dr {params.dist_range} -tc {threads} '
        '-o {output.hla_bam_R1} {params.hla_ref} {input.fastq_R1} ; '
        '/usr/local/bin/razers3 -i {params.perc_indent} -m {params.max_hits} -dr {params.dist_range} -tc {threads} '
        '-o {output.hla_bam_R2} {params.hla_ref} {input.fastq_R2} '

# optitype samtools conversion step
rule optitype_convert_extracted:
    input:
        hla_bam_R1='main_run/{patient}/{sample}/outputs/optitype_razers/sample_R1.bam',
        hla_bam_R2='main_run/{patient}/{sample}/outputs/optitype_razers/sample_R2.bam',
    output:
        hla_fastq_R1='main_run/{patient}/{sample}/outputs/optitype_convert_extracted/sample_R1.fastq',
        hla_fastq_R2='main_run/{patient}/{sample}/outputs/optitype_convert_extracted/sample_R2.fastq',
    singularity:
        "docker://biocontainers/samtools:v1.9-4-deb_cv1",
    shell:
        # create the command line using the input
        'samtools bam2fq {input.hla_bam_R1} > {output.hla_fastq_R1}; '
        'samtools bam2fq {input.hla_bam_R2} > {output.hla_fastq_R2}'

# optitype step
rule optitype:
    input:
        hla_fastq_R1='main_run/{patient}/{sample}/outputs/optitype_convert_extracted/sample_R1.fastq',
        hla_fastq_R2='main_run/{patient}/{sample}/outputs/optitype_convert_extracted/sample_R2.fastq',
    output: # main_run/A002/ADT002/outputs/optitype/sample/sample_coverage_plot.pdf
        coverage='main_run/{patient}/{sample}/outputs/optitype/sample/sample_coverage_plot.pdf',
        hla='main_run/{patient}/{sample}/outputs/optitype/sample/sample_result.tsv',
    singularity:
        "docker://fred2/optitype",
    shell:
        'OPTITYPE_DIR=/usr/local/bin/OptiType; ' 
        # create the command line using the input paths
        'python2.7 $OPTITYPE_DIR/OptiTypePipeline.py -i {input.hla_fastq_R1} {input.hla_fastq_R2} ' 
        # add the output path and some parameters
        '--dna --config $OPTITYPE_DIR/config.ini --prefix sample -v -o {params.outdir}' 

rule maf_to_vcf:
    input:
        maf=_get_maf, 
        ref=config['reference_fasta'],
    output:
        vcf='main_run/{patient}/{sample}/outputs/pvacseq_input/consensus.vcf',
    singularity:
        '/juno/work/shah/mondrian/singularity/variant_v0.0.26.sif', # contains maf2vcf.pl
    shell:
        'maf2vcf.pl --input-maf {input.maf} --ref-fasta {input.ref} ' 
        '--output-dir {params.outdir} ' 
        '--output-vcf {output}'

rule vep:
    input:
        vcf='main_run/{patient}/{sample}/outputs/pvacseq_input/consensus.vcf',
        ref=config['vep']['fasta'],
    output:
        vcf='main_run/{patient}/{sample}/outputs/pvacseq_input/consensus.vep.vcf',
    params:
        vep_run_file=config['vep']['run_file'],
        vep_plugins_dir=config['vep']['plugins_dir'],
        vep_cache_dir=config['vep']['cache_dir'],
    shell:
        'module load perl && ' 
        '{params.vep_run_file} '
        '--input_file {input.vcf} --format vcf '
        '--output_file {output.vcf} --vcf '
        '--symbol --terms SO --tsl --hgvs ' 
        '--fasta {input.ref} ' 
        '--offline --cache ' 
        '--plugin Frameshift --plugin Wildtype ' 
        '--dir_plugins {params.plugins_dir} '
        '--dir_cache {params.cache_dir} '
        '--force_overwrite'

rule vt_decompose:
    input:
        vcf='main_run/{patient}/{sample}/outputs/pvacseq_input/consensus.vep.vcf',
    output:
        vcf='main_run/{patient}/{sample}/outputs/pvacseq_input/consensus.vep.decomposed.vcf',
    singularity:
        "docker://zlskidmore/vt",
    shell:
        'vt decompose '
        '-s {input.vcf} '
        '-o {output.vcf}' 

rule bam_readcount:
    input:
        vcf='main_run/{patient}/{sample}/outputs/pvacseq_input/consensus.vep.decomposed.vcf',
        bam=_get_bam,
        ref=config['reference_fasta'],
    output:
        snv_tsv='main_run/{patient}/{sample}/outputs/pvacseq_input/{sample}_bam_readcount_snv.tsv',
        indel_tsv='main_run/{patient}/{sample}/outputs/pvacseq_input/{sample}_bam_readcount_indel.tsv',
    params:
        outdir='main_run/{patient}/{sample}/outputs/pvacseq_input'
    singularity:
        '/juno/work/shah/users/chois7/apollo/neoantigen/bam-readcount_helper_v0.0.3.sif',
    shell:
        'python /usr/bin/bam_readcount_helper.py ' 
        '{input.vcf} ' 
        '{wildcards.sample} '
        '{input.ref} ' 
        '{input.bam} '
        '{params.outdir}'
    
rule vcf_readcount_annotator_snv:
    input:
        vcf='main_run/{patient}/{sample}/outputs/pvacseq_input/consensus.vep.decomposed.vcf',
        snv_tsv='main_run/{patient}/{sample}/outputs/pvacseq_input/{sample}_bam_readcount_snv.tsv',
    output:
        vcf='main_run/{patient}/{sample}/outputs/pvacseq_input/consensus.vep.decomposed.snv_annotated.vcf',
    singularity:
        'docker://griffithlab/vatools',
    shell:
        'vcf-readcount-annotator '
        '{input.vcf} {input.snv_tsv} '
        'DNA -s {wildcards.sample} '
        '-t snv -o {output.vcf}'

rule vcf_readcount_annotator_indel:
    input:
        vcf='main_run/{patient}/{sample}/outputs/pvacseq_input/consensus.vep.decomposed.snv_annotated.vcf',
        indel_tsv='main_run/{patient}/{sample}/outputs/pvacseq_input/{sample}_bam_readcount_indel.tsv',
    output:
        vcf='main_run/{patient}/{sample}/outputs/pvacseq_input/consensus.vep.decomposed.both_annotated.vcf',
    singularity:
        'docker://griffithlab/vatools',
    shell:
        'vcf-readcount-annotator '
        '{input.vcf} {input.indel_tsv} '
        'DNA -s {wildcards.sample} '
        '-t indel -o {output.vcf}'

rule pvacseq:
    input:
        vcf='main_run/{patient}/{sample}/outputs/pvacseq_input/consensus.vep.decomposed.both_annotated.vcf',
        hla='main_run/{patient}/{sample}/outputs/optitype/sample/sample_result.tsv',
    output: # main_run/A002/ADT002/outputs/pvacseq/MHC_Class_I/ADT002.filtered.tsv
        tsv='main_run/{patient}/{sample}/outputs/pvacseq/MHC_Class_I/{sample}.filtered.tsv',
    params:
        outdir='main_run/{patient}/{sample}/outputs/pvacseq',
        # TODO: reduce modes; e.g. use NetMHCcons for NetMHC, (NetMHCpan), and PickPocket
        #modes = 'MHCflurry MHCnuggetsI MHCnuggetsII NNalign NetMHC PickPocket SMM SMMPMBEC SMMalign',
        modes = 'MHCflurry MHCnuggetsI MHCnuggetsII NNalign NetMHC PickPocket SMM SMMPMBEC SMMalign',
    singularity: '/juno/work/shah/users/chois7/apollo/pvactools_latest.sif',
    resources:
        mem_mb=config['pvacseq']['mem_mb'],
        disk_mb=config['pvacseq']['disk_mb'],
    shell:
        ## input.hla content:
        #         A1      A2      B1      B2      C1      C2      Reads   Objective
        # 0       A*24:02 A*32:01 B*14:01 B*15:01 C*03:03 C*08:02 1300.0  1241.4799999999998
        'hlas=$(cat {input.hla} | tail -n 1 | cut -f2,3,4,5,6,7) && '
        'hlacs="" && '
        'for hla in ${{hlas}}; do echo $hla; hlacs+="HLA-${{hla}},"; done && '
        'pvacseq run ' 
        '--iedb-install-directory /opt/iedb ' # local IEDB directory 
        '--blastp-path /opt/ncbi-blast-2.12.0+/bin/blastp ' # local BLASTP binary #'--keep-tmp-files ' + # keep temp
        '-t 8 ' # threading
        '{input.vcf} '
        '{wildcards.sample} '
        '${{hlacs::-1}} ' # HLA-A*24:02,...,HLA-C*08:02, -> rm comma at the end
        '{params.modes} '  
        '{params.outdir}'

rule pvacseq_expression_filter: 
    input:
        tsv='main_run/{patient}/{sample}/outputs/pvacseq/MHC_Class_I/{sample}.filtered.tsv',
        kallisto=_get_kallisto,
    output:
        tsv='main_run/{patient}/{sample}/outputs/pvacseq/MHC_Class_I/{sample}.expression_filtered.tsv',
    params:
        cutoff = 1.0,
    singularity: 'docker://soymintc/clickpdvcf',
    shell:
        'python scripts/expression_filter.py '
        '--pvacseq {input.tsv} '
        '--kallisto {input.kallisto} '
        '--cutoff {params.cutoff} '
        '--output {output.tsv}'

