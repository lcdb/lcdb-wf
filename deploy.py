#!/usr/bin/env python

import os
import sys

import tempfile
import argparse
import subprocess as sp
import datetime
import json
import fnmatch
import logging
import hashlib
from pathlib import Path
from distutils import filelist

# Determine default staging area, used in help
default_staging = "/tmp/{0}-lcdb-wf-staging".format(os.getenv('USER'))

usage = f"""
This script assists in the deployment of relevant code from the lcdb-wf
repository to a new deployment directory for running an analysis. It is
intended to be run in a standalone fashion such that with just the script you
can download and deploy a specified version of the workflows.

For example, the following command will clone the GitHub repo to {default_staging},
check out the v9.999 branch, copy the files needed for RNA-seq over to the
"my_analysis_dir" directory, store a read-only file .lcdb-wf-deployment.yaml
with the metadata of the repo used for cloning, and build the conda
environments within "my_analysis_dir":

    ./deploy.py \\
        --clone \\
        --dest my_analysis_dir \\
        --flavor rnaseq \\
        --build-envs \\
        --branch v9.999

Compared to directly cloning the repo, this results in a cleaner deployment
directory that does not have various test infrastructure or workflows not
relevant to the project.
"""

logging.basicConfig(
    format="%(asctime)s [%(module)s] %(message)s",
    level=logging.DEBUG,
    datefmt="%Y-%m-%dT%H:%M:%S",
)

# ANSI color escape codes
DEFAULT = "\x1b[39m"
RED = "\x1b[31m"
GREEN = "\x1b[32m"
YELLOW = "\x1b[33m"
GRAY = "\x1b[37m"
WHITE = "\x1b[97m"
BLUE = "\x1b[34m"
RESET = "\x1b[0m"


def debug(s):
    logging.debug(GRAY + s + RESET)


def info(s):
    logging.info(GREEN + s + RESET)


def warning(s):
    logging.warning(YELLOW + s + RESET)


def error(s):
    logging.error(RED + s + RESET)


def write_include_file(source, flavor='all'):

    # Patterns follow that of MANIFEST.in
    # (https://packaging.python.org/en/latest/guides/using-manifest-in/),
    # and distutils.filelist is used below to parse them.

    PATTERN_DICT = {
        'rnaseq': [
            'include workflows/rnaseq/Snakefile',
            'recursive-include workflows/rnaseq/config *',
            'include workflows/rnaseq/rnaseq_trackhub.py',
            'recursive-include workflows/rnaseq/downstream *.Rmd',
            'recursive-include workflows/rnaseq/downstream *.yaml',
        ],
        'chipseq': [
            'include workflows/chipseq/Snakefile',
            'recursive-include workflows/chipseq/config *',
            'include workflows/chipseq/chipseq_trackhub.py',
        ],
        'all': [
            'recursive-include wrappers *',
            'recursive-include include *',
            'recursive-include lib *', 
            'include env.yml env-r.yml .gitignore',
            'include workflows/references/Snakefile',
            'recursive-include workflows/references/config *',
            'global-exclude __pycache__',
        ],
        'full': [
            'include workflows/colocalization/Snakefile',
            'recursive-include workflows/colocalization/config *',
            'recursive-include workflows/colocalization/scripts *',
            'recursive-include workflows/figures *',
            'recursive-include workflows/external *',
        ]

    }

    patterns = []
    if flavor in ('full', 'rnaseq'):
        patterns.extend(PATTERN_DICT['rnaseq'])
    if flavor in ('full', 'chipseq'):
        patterns.extend(PATTERN_DICT['chipseq'])
    if flavor == 'full':
        patterns.extend(PATTERN_DICT['full'])
    patterns.extend(PATTERN_DICT['all'])

    def fastwalk(path):
        """
        Find all files recursively, but short-circuit if we get to a conda env to
        avoid traversing those many files.
        """
        path = str(path)
        for root, dirs, files in os.walk(path, topdown=True):
            if 'conda-meta' in dirs:
                dirs[:] = []
                files[:] = []
            for d in dirs:
                yield os.path.join(root, d).replace(path + '/', '')
            for f in files:
                yield os.path.join(root, f).replace(path + '/', '')

    f = filelist.FileList()
    f.allfiles = list(fastwalk(source))
    for pattern in patterns:
        f.process_template_line(pattern)
    f.sort()
    f.remove_duplicates()

    under_version_control = sorted(
        sp.check_output(
            ["git", "ls-tree", "-r", "HEAD", "--name-only"],
            universal_newlines=True,
            cwd=source,
        ).splitlines(False),
    )

    to_transfer = list(set(under_version_control).intersection(f.files))
    include = tempfile.NamedTemporaryFile(delete=False).name
    with open(include, 'w') as fout:
        fout.write('\n\n')
        fout.write('\n'.join(to_transfer))

    return include


def clone_repo(dest, branch="master", mismatch_ok=False):

    if Path(dest).exists():
        error("Path {dest} already exists, aborting!".format(**locals()))
        sys.exit(1)

    URL = "https://github.com/lcdb/lcdb-wf.git"
    info("cloning {URL} to {dest}".format(**locals()))
    cmds = ["git", "clone", URL, dest]
    sp.check_call(cmds)
    sp.check_call(["git", "checkout", branch], cwd=dest)
    info("Cloned to {dest} and will deploy using branch '{branch}'".format(**locals()))

    # check to see if this very file that is running is the same as the one
    # that was just cloned -- otherwise it's out of date.

    def check_md5(f):
        contents = open(f).read()
        contents = str(contents).encode("utf-8")
        return hashlib.md5(contents).hexdigest()

    this_md5 = check_md5(__file__)
    that_md5 = check_md5(os.path.join(dest, "deploy.py"))

    if (this_md5 != that_md5) and not mismatch_ok:
        full_here = Path(__file__).resolve()
        full_there = Path(dest) / "deploy.py"
        error(
            "Files {full_here} and {full_there} do not match! ".format(**locals()) +
            "The deploy script you are running appears to be out of date. "
            "Please get an updated copy from https://github.com/lcdb/lcdb-wf, perhaps "
            "with 'wget https://raw.githubusercontent.com/lcdb/lcdb-wf/master/deploy.py'"
        )
        sys.exit(1)


def rsync(include, source, dest, rsync_args):
    rsync = [
        "rsync",
        "--relative",
        rsync_args,
        "--files-from={}".format(include),
        source,
        dest,
    ]
    sp.check_call(rsync)


def deployment_json(source, dest):
    """
    Build the .lcdb-wf-deployment.json file
    """
    # First, the last commit:
    commit, message = (
        sp.check_output(
            ["git", "log", "--oneline", "-1"], universal_newlines=True, cwd=source
        )
        .strip()
        .split(" ", 1)
    )

    # When we're deploying:
    now = datetime.datetime.strftime(datetime.datetime.now(), "%Y%m%d%H%M")

    # Where the remote was:
    remotes = sp.run(
        ["git", "remote", "-v"],
        universal_newlines=True,
        cwd=source,
        check=True,
        stdout=sp.PIPE,
        stderr=sp.STDOUT,
    )
    remotes = [i.strip() for i in remotes.stdout.splitlines()]

    # The branch we're deploying from:
    branch = sp.run(
        ["git", "branch"],
        universal_newlines=True,
        cwd=source,
        check=True,
        stdout=sp.PIPE,
        stderr=sp.STDOUT,
    )
    branch = [i for i in branch.stdout.splitlines() if i.startswith("*")]
    assert len(branch) == 1
    branch = branch[0]
    branch = branch.split("* ")[1]

    d = {
        "git": {
            "commit": commit,
            "message": message,
            "remotes": remotes,
            "branch": branch,
        },
        "timestamp": now,
    }
    log = Path(dest) / ".lcdb-wf-deployment.json"

    with open(log, "w") as fout:
        fout.write(json.dumps(d) + "\n")
    os.chmod(log, 0o440)

    info("Wrote details of deployment to {log}".format(**locals()))


def build_envs(dest, additional_main=None, additional_r=None, conda_frontend="conda"):
    """
    Build conda environments.

    Parameters
    ----------

    dest : str
        Destination path. This is the project directory (likely as specified on
        the command line with --dest) in which the env and env-r yaml files
        should already exist. Envs will be created in here.

    additional_main : list
        Other packages to install, e.g., a snakemake plugin needed for
        a cluster profile, into the main environment.

    additional_r : list
        Other packages to install into the R environment.

    conda_frontend : 'mamba' | 'conda'
        Which front-end to use (terminology borrowed from Snakemake)

    """
    mapping = [
        ("./env", "env.yml", additional_main),
        ("./env-r", "env-r.yml", additional_r),
    ]
    for env, yml, additional in mapping:
        info("Building environment " + os.path.join(dest, env))
        if additional:
            info(f"Adding {additional} to environment")

        try:
            # conda and mamba can be hard to kill, possibly because they're
            # doing threaded things. So we use Popen explicitly to capture the
            # process ID so it can be killed if the user hits ^C.
            #
            # Here we use Popen directly in order to get the process ID so that
            # it can be explictly kill upon KeyboardInterrupt.
            cmds = [
                conda_frontend,
                "env",
                "create",
                "-p",
                env,
                "--file",
                yml,
            ]
            if additional:
                cmds += additional
            p = sp.Popen(cmds, universal_newlines=True, cwd=dest)
            p.wait()

        except KeyboardInterrupt:
            print("")
            error("Killing running {conda_frontend} job, '".format(**locals()) + " ".join(cmds))
            p.kill()
            sys.exit(1)

        if p.returncode:
            error("Error running {conda_frontend}, '".format(**locals()) + " ".join(cmds))
            sys.exit(1)

        full_env = Path(dest) / env
        info("Created env {full_env}".format(**locals()))


if __name__ == "__main__":

    ap = argparse.ArgumentParser(usage=usage)
    ap.add_argument(
        "--flavor",
        default="full",
        help="""Options are {0}. Default is full.""".format(['full', 'rnaseq', 'chipseq']),
    )
    ap.add_argument(
        "--dest", help="""Destination directory in which to copy files""", required=True
    )

    ap.add_argument(
        "--clone",
        action="store_true",
        help=f"""Make a new clone to a staging area (at the location specified
        by --staging which defaults to {default_staging}) and deploy from
        there. Useful if using this script as a standalone tool. You can also
        use --branch to configure which branch to deploy from that clone."""
    )

    ap.add_argument(
        "--staging",
        help="""Only used when --clone is specified. Clone the main git repo to
        this directory and do a diff on the deploy.py script found there to
        ensure this one is up-to-date, and if so then proceed using the new clone as the source.
        """,
    )

    ap.add_argument(
        "--branch",
        help="Branch to checkout if using --staging to clone a temporary copy. Default is %(default)s.",
        default="master",
    )

    ap.add_argument(
        "--build-envs",
        action="store_true",
        help="""If specified, conda environments with all dependencies will be
        installed into directories called "env" and "env-r" within the directory
        provided for --dest.""",
    )
    ap.add_argument(
        "--conda-frontend",
        help="Set program (conda or mamba) to use when creating environments. Default is %(default)s.",
        default="conda",
    )
    ap.add_argument(
        "--rsync-args",
        help="Options for rsync when deploying to a new directory. Default is %(default)s.",
        default="-rlt"
    )

    ap.add_argument(
        "--additional-main",
        help="""Additional packages to install in main environment (only
        relevant with --build-envs). For example,
        'snakemake-executor-plugin-cluster-generic' to support a cluster
        profile.""",
        nargs="+"
    )
    ap.add_argument(
        "--additional-r",
        help="Additional packages to install in R environment (only relevant with --build-envs)",
        nargs="+"
    )

    ap.add_argument(
        "--mismatch-ok",
        action="store_true",
        help="Used for testing")
    args = ap.parse_args()
    dest = args.dest
    flavor = args.flavor

    if args.staging and not args.clone:
            print("ERROR: --staging was specified but --clone was not. Did you want to use --clone?", file=sys.stderr)
            sys.exit(1)
    if args.clone:
        if args.staging is None:
            args.staging = default_staging
        source = os.path.abspath(args.staging)
        clone_repo(args.staging, args.branch, mismatch_ok=args.mismatch_ok)
    else:
        source = Path(__file__).parent.resolve()

    include = write_include_file(source, flavor)
    rsync(include, source, dest, args.rsync_args)
    deployment_json(source, dest)

    if args.additional_main and additional_main_from_env_var:
        print(
            "ERROR: Unset LCDBWF_ADDITIONAL_MAIN env var if you want to use the --additional-main argument."
        )
        sys.exit(1)

    if additional_main_from_env_var:
        if args.additional_main:
            print(
                "ERROR: Unset LCDBWF_ADDITIONAL_MAIN env var if you want to use the --additional-main argument."
            )
            sys.exit(1)
        additional_main = [additional_main_from_env_var]
    else:
        additional_main = args.additional_main

    if args.build_envs:
        build_envs(
            dest,
            additional_main=additional_main,
            additional_r=args.additional_r,
            conda_frontend=args.conda_frontend,
        )

    warning("Deployment complete in {args.dest}".format(**locals()))
