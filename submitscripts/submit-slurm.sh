#!/bin/bash

snakemake   \
        --jobs 100 -s Snakefile.smk --use-conda --cluster-config submitscripts/cluster.yaml --cluster "sbatch --parsable --partition={cluster.queue} --job-name=asv.{rule}.{wildcards} --mem={cluster.mem}gb --time={cluster.time} --ntasks={cluster.threads} --nodes={cluster.nodes} --verbose"

