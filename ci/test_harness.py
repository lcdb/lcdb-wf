#!/usr/bin/env python

import argparse
from snakemake import shell
from pathlib import Path
import json

ap = argparse.ArgumentParser()
ap.add_argument('name', help='Name of test to run')
ap.add_argument('--list-names', action='store_true', help='List names of available tests and then exit')
ap.add_argument('--deploy', help='Location in which to deploy')
ap.add_argument('--env', help='Environment to source before running')
ap.add_argument('--threads', help='Number of threads to pass to snakemake, default is %(default)s', default=1)
args = ap.parse_args()

REPO = Path(__file__).parent.resolve()

env = Path(args.env).resolve()

info = json.loads(shell.check_output('conda info --json'))
conda_prefix = info['conda_prefix']
shell(
    'set +eu; source {conda_prefix}/bin/activate {env}; set -eu; '
    'cd {args.deploy}; '
    'cp -r workflows/rnaseq workflows/rnaseq-misc-test; '
    'cp -r workflows/rnaseq/data /tmp/data; '
    'cd workflows/rnaseq-misc-test; '
    'rm -r data; '
    './run_test -j {args.threads} --use-conda -k -p -r --until cutadapt '
    '--configfile {REPO}/test/test_configs/test_rnaseq_config.yaml '
    '--config sampletable={REPO}/test/test_configs/test_sra_sampletable.tsv '
)
