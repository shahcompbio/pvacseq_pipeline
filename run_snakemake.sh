CLUSTER_CMD=("bsub -n {threads} -R {cluster.resources} -M {cluster.memory} -o {cluster.output} -J {cluster.name} -W {cluster.time}") 
snakemake --configfile config.yaml \
    --use-singularity \
    --singularity-args "--bind /juno --bind /rtsess01" \
    --jobs 30 \
    --skip-script-cleanup \
    --cluster-config ./cluster.yaml \
    --cluster "${CLUSTER_CMD}" \
    --cluster-cancel bkill \
    --rerun-incomplete #--dry-run
