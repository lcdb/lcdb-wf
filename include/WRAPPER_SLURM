#!/bin/bash
#SBATCH --job-name="lcdb-wf"
#SBATCH --partition="norm"
#SBATCH --time=24:00:00

# make logdir
if [[ ! -e logs ]]; then mkdir -p logs; fi

# Run snakemake
(
    time snakemake \
    -p \
    --directory $PWD \
    -k \
    --rerun-incomplete \
    --jobname "s.{rulename}.{jobid}.sh" \
    -j 999 \
    --cluster-config config/clusterconfig.yaml \
    --verbose \
    --cluster 'sbatch {cluster.prefix} --cpus-per-task={threads}  --output=logs/{rule}.o.%j --error=logs/{rule}.e.%j' \
    --use-conda \
    --configfile config/config.yaml \
    --latency-wait=300 \
    ) > "Snakefile.log" 2>&1

SNAKE_PID=$!

finish(){
    echo 'Stopping running snakemake job.'
    kill -SIGINT $SNAKE_PID
    exit 0
}
trap finish SIGTERM

wait $SNAKE_PID
