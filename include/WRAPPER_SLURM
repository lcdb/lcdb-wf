#!/bin/bash
#SBATCH --job-name="lcdb-wf"
#SBATCH --partition="norm"
#SBATCH --time=12:00:00
#SBATCH --gres=lscratch:5

# make logdir
if [[ ! -e logs ]]; then mkdir -p logs; fi

# use SNAKEMAKE_PROFILE if one is set in the environment variables
if [ -z "$LCDBWF_SNAKEMAKE_PROFILE" ]; then
    if [ -z "$SNAKEMAKE_PROFILE" ]; then
        # no snakemake profile found
        PROFILE_CMD=""
        echo "No environment variable SNAKEMAKE_PROFILE or LCDBWF_SNAKE_PROFILE found."
        echo "snakemake will run in single job."
    else
        # generic SNAKEMAKE_PROFILE found
        PROFILE_CMD="--profile $SNAKEMAKE_PROFILE"
    fi
else
# LCDBWF_SNAKEMAKE_PROFILE found, this takes priority if both profile variables are set
PROFILE_CMD="--profile $LCDBWF_SNAKEMAKE_PROFILE"
fi

# Run snakemake
(
    time snakemake \
    -p \
    --directory $PWD \
    -k \
    --restart-times 3 \
    --rerun-incomplete \
    --jobname "s.{rulename}.{jobid}.sh" \
    -j 999 \
    --use-conda \
    --configfile config/config.yaml \
    $PROFILE_CMD \
    "$@"
    ) > "Snakefile.log" 2>&1

SNAKE_PID=$!

finish(){
    echo 'Stopping running snakemake job.'
    kill -SIGINT $SNAKE_PID
    exit 0
}
trap finish SIGTERM

wait $SNAKE_PID
