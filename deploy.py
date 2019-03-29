import os
import tempfile
import argparse
import subprocess as sp
import datetime
import json

HERE = os.path.dirname(__file__)

usage = """
This script assists in the deployment of lcdb-wf to working directories.

The lcdb-wf repository contains infrastructure for testing that is not
typically needed when using it in practice. Furthermore, you might not need all
possible workflows.

This script copies over only the files requred for each "flavor" of analysis
(rnaseq, chipseq, colocalization, full) and also stores a file,
`.lcdb-wf-deployment.yaml`, containing details about the git commit that was
used and the timestamp. This can be used to compare changes and stay
up-to-date.
"""

ap = argparse.ArgumentParser(usage=usage)
ap.add_argument('--flavor', default='full', help='''Options are rnaseq, chipseq, colocalization, full. Default is full.''')
ap.add_argument('--dest', help='''Destination directory in which to copy files''')
args = ap.parse_args()
dest = args.dest
flavor = args.flavor

flavors = {
    'all': {
        'include': [
            'wrappers/wrappers',
            'include',
            'lib',
            'requirements.txt',
        ],
        'exclude': [
            'wrappers/wrappers/demo',
            'workflows/*/run_test.sh',

            # The following files to exclude are those that are created from
            # a test run.
            'lib/__pycache__',
            'lib/postprocess/__pycache__',
            'include/AnnotationHubCache',
            'workflows/*/Snakefile.test',
            'workflows/*/references_data',
            'workflows/*/.snakemake',
            'workflows/*/data',
            'workflows/rnaseq/downstream/rnaseq_cache',
            'workflows/rnaseq/downstream/rnaseq_files',
            'workflows/rnaseq/downstream/final_clusters',
            'workflows/rnaseq/downstream/*.tsv*',
            'workflows/rnaseq/downstream/*log',
            'workflows/rnaseq/downstream/*html',
            'workflows/colocalization/results',
        ],
    },
    'chipseq': [
        'workflows/chipseq',
        'workflows/references',
    ],
    'rnaseq': [
        'workflows/rnaseq',
        'workflows/rnaseq/downstream/',
        'workflows/references',
    ],
    'colocalization': [
        'workflows/colocalization',
    ],

    'full': [
        'workflows',
    ],
}

paths = set(flavors['all']['include'])
paths = paths | set(flavors[flavor])

exclude = tempfile.NamedTemporaryFile(delete=False).name
with open(exclude, 'w') as fout:
    fout.write('\n'.join(flavors['all']['exclude']))

include = tempfile.NamedTemporaryFile(delete=False).name
with open(include, 'w') as fout:
    fout.write('\n'.join(paths) + '\n')


sp.check_call([
    'rsync',
    '--relative',
    '-ar',
    '--files-from={}'.format(include),
    '--exclude-from={}'.format(exclude),
    HERE,
    dest])


commit, message = sp.check_output(
    ['git', 'log', '--oneline', '-1'],
    universal_newlines=True
).strip().split(' ', 1)
now = datetime.datetime.strftime(datetime.datetime.now(), '%Y%m%d%H%M')
remotes = sp.check_output(
    ['git', 'remote', '-v'],
    universal_newlines=True
)
remotes = [i.strip() for i in remotes.splitlines()]
branch = sp.check_output([
    'git', 'branch'], universal_newlines=True)
branch = [i for i in branch.splitlines() if i.startswith('*')]
assert len(branch) == 1
branch = branch[0]
branch = branch.split('* ')[1]

d = {
    'git': {
        'commit': commit,
        'message': message,
        'remotes': remotes,
        'branch': branch,
    },
    'timestamp': now}
log = os.path.join(dest, '.lcdb-wf-deployment.json')
with open(log, 'w') as fout:
    fout.write(json.dumps(d) + '\n')
os.chmod(log, 0o440)
