#!/bin/bash
module load snakemake/5.19.2 python/3.7
snakemake -s annotation_pipeline.snake --jobs 200 --cluster '{cluster.scheduler} {cluster.queue} {cluster.export_env} {cluster.cwd} {cluster.mem} {cluster.n_cpu}{cluster.threads}' --cluster-config Additional_files/config_cluster_slurm.yaml --use-singularity $1 --singularity-args '--home ~/' | tee stdout.txt
