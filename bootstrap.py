#!/usr/bin/env python
import os
import sys
import argparse
import tempfile
from snakemake.shell import shell
import glob

usage = """
Clones the repo into a new directory ready for a new analysis.

Specifically:

    - clone current repo to a temp dir
    - check out a new branch in the tmp repo
    - git rm everything that isn't listed on manifest.txt and commit on the new branch
    - add the github repo as "upstream" remote
    - copy the directory over, but warn about conflicts
"""
BLUE = '\033[94m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
RED = '\033[91m'
END = '\033[0m'

ap = argparse.ArgumentParser()
ap.add_argument('dest', help='Target directory. Will be created if needed')
ap.add_argument('--label', help='creates a new branch of this name.', required=True)
ap.add_argument('--https', action='store_true',
                help='''when setting the upstream remote, use https rather than
                'ssh url''')
args = ap.parse_args()

source_dir = os.path.dirname(os.path.realpath(__file__))
dest_dir = args.dest

if not os.path.exists(dest_dir):
    os.makedirs(dest_dir)
else:
    if os.path.exists(os.path.join(dest_dir, '.git')):
        print(RED + "\nERROR: {dest_dir} appears to be an existing git repo, aborting.\n".format(**locals()) + END)
        sys.exit(1)


def recursive_find(paths):
    files = []
    for path in paths:
        if not os.path.isdir(path):
            files.append(path)
            continue
        for root, dirnames, filenames in os.walk(path):
            for filename in filenames:
                files.append(os.path.join(root, filename))
    return files


def _run(cmd):
    print('\t' + YELLOW + cmd + END)
    shell(cmd)

tmpdir = tempfile.mkdtemp()
_run(
    'cd {tmpdir} '
    '&& git clone {source_dir} . '
    '&& git submodule init '
    '&& git submodule update '
    '&& git checkout -b {args.label}'.format(**locals()))




# Whitelist of globs indicating files to keep
MANIFEST = [i.strip() for i in open('manifest.txt') if
            len(i.strip()) > 0 and not i.startswith('#')]

keep = []
for m in MANIFEST:
    keep.extend(recursive_find(glob.glob(os.path.join(tmpdir, m))))

# exhaustive list of files in the repo.
existing_in_repo = recursive_find([tmpdir])

# Anything not in the whitelist gets deleted from the temp repo
to_remove = set(existing_in_repo).difference(set(keep))

# If we have something existing in the dest (say, an existing config file),
# then add that to the remove list so that when copying over we don't
# overwrite.
existing_in_dest = recursive_find([args.dest])
to_remove = to_remove.union(existing_in_dest)

# add the remote
if not args.https:
    remote = 'git@github.com:lcdb/lcdb-wf.git'
else:
    remote = 'https://github.com/lcdb/lcdb-wf.git'
_run('cd {tmpdir} && git remote add upstream {remote}'.format(**locals()))


# Now remove anything specified 
for r in list(to_remove):
    r = os.path.relpath(r, tmpdir)
    _run('(cd {tmpdir} && git rm -r {r} --quiet)'.format(**locals()))

_run('cd {tmpdir} && git commit -a -m "cleanup for new project \'{args.label}\' in {dest_dir}" --quiet'.format(**locals()))

_run('rsync -r {tmpdir}/ {dest_dir}/'.format(**locals()))

