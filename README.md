# Bulk RNA-seq variant discovery pipeline
Follows pVACseq DNA-Seq pipeline as in [link](https://pvactools.readthedocs.io/en/latest/pvacseq.html)

## Prerequisites
- You need `snakemake` to be installed
- Data will be stored in `./main_run`, but you can change the `Snakefile` and `rules/*.smk` accordingly to change your results path.

## Running on APOLLO cohort
```bash
bash run_snakemake.sh
```
