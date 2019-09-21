import os
import tempfile
import argparse
import subprocess as sp
import datetime
import json
import fnmatch


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
ap.add_argument('--flavor', default='full', help='''Options are rnaseq,
                chipseq, colocalization, full. Default is full.''')
ap.add_argument('--dest', help='''Destination directory in which to copy files''')
ap.add_argument('--build-env', action='store_true', help='''If specified,
a conda environment with all dependencies will be installed into a directory
                called "env" within the directory provided for --dest.  ''')

args = ap.parse_args()
dest = args.dest
flavor = args.flavor


# Inclusions and exclusions:
# Excluding files by name or by location is easy: --exclude=*~,
# --exclude=/some/relative/location (relative to the source argument, e.g. this
# excludes ~/LaTeX/some/relative/location).
# If you only want to match a few files or locations, include them, include
# every directory leading to them (for example with --include=*/), then exclude
# the rest with --exclude='*'. This is because:
# If you exclude a directory, this excludes everything below it. The excluded
# files won't be considered at all.
# If you include a directory, this doesn't automatically include its contents.
# In recent versions, --include='directory/***' will do that.
# For each file, the first matching rule applies (and anything never matched is
# included).
# Patterns:
# If a pattern doesn't contain a /, it applies to the file name sans directory.
# If a pattern ends with /, it applies to directories only.
# If a pattern starts with /, it applies to the whole path from the directory
# that was passed as an argument to rsync.
# * any substring of a single directory component (i.e. never matches /); **
# matches any path substring.
flavors = {
    'all': {
        'include': [
            'wrappers/wrappers',
            'include',
            'lib',
            'requirements.txt',
            '.gitignore',
        ],
        'exclude': [

            '/.buildkite/*',
            '/ci/*',
            '/.circleci/*',
            '/config/*',
            '/deploy.py',
            '/docs/***',
            '/include/AnnotationHubCache',
            '/lib/postprocess/__pycache__',
            '/lib/__pycache__',
            '*/.pytest_cache/*',
            '/README.md',
            '/test/*',
            '/.travis.yml',
            '/workflows/*/results',
            '/workflows/*/data',
            '/workflows/figures/*',
            '/workflows/*/references_data',
            '/workflows/*/references_dir',
            '/workflows/*/reports',
            '/workflows/rnaseq/downstream/final_clusters',
            '/workflows/rnaseq/downstream/*html',
            '/workflows/rnaseq/downstream/*log',
            '/workflows/rnaseq/downstream/rnaseq_cache',
            '/workflows/rnaseq/downstream/rnaseq_files',
            '/workflows/rnaseq/downstream/*.tsv*',
            '/workflows/*/run_test.sh',
            '/workflows/*/Snakefile.test',
            '/workflows/*/.snakemake',
            '/wrappers/demo/*',
            '/wrappers/test/*',
            '/wrappers/test_toy.py',

        ],
    },
    'chipseq': [
        'workflows/chipseq/*',
        'workflows/references/*',
    ],
    'rnaseq': [
        'workflows/rnaseq/*',
        'workflows/rnaseq/config/*',
        'workflows/rnaseq/downstream/*',
        'workflows/references/*',
    ],
    'colocalization': [
        'workflows/colocalization/*',
    ],

    'full': [
        'workflows/*',
    ],
}


under_version_control = sorted(sp.check_output(
    [ 'git', 'ls-tree', '-r', 'HEAD', '--name-only'],
    universal_newlines=True).splitlines(False))


def filter_out_excluded(filenames, patterns_to_exclude):
    keep = []
    for filename in filenames:
        if not any(fnmatch.fnmatch(filename, pattern) for pattern in patterns_to_exclude):
            keep.append(filename)
    return keep

def filter_out_other_workflows(filenames, patterns_to_include, prefix='workflows'):
    keep = []
    for filename in filenames:
        if not filename.startswith(prefix):
            keep.append(filename)
            continue
        if any(fnmatch.fnmatch(filename, pattern) for pattern in patterns_to_include):
            keep.append(filename)
    return keep

keep = filter_out_excluded(under_version_control, flavors['all']['exclude'])

keep = filter_out_other_workflows(keep, flavors[flavor])

exclude = tempfile.NamedTemporaryFile(delete=False).name
exclude = '.exclude'
with open(exclude, 'w') as fout:
    fout.write('\n'.join(flavors['all']['exclude']))

include = tempfile.NamedTemporaryFile(delete=False).name
include = '.include'
with open(include, 'w') as fout:
    fout.write('\n\n')
    fout.write('\n'.join(keep) + '\n')

sp.check_call([
    'rsync',
    '--relative',
    '-v',
    '-ar',
    '--progress',
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

if args.build_env:
    sp.check_call(
        ['conda', 'create', '-y', '-p', './env', '--file', 'requirements.txt', '-c',
         'conda-forge', '-c', 'bioconda', '-c', 'defaults'],
        universal_newlines=True,
        cwd=dest)
