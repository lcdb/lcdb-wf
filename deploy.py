#!/usr/bin/env python

import os
import sys

try:
    from pathlib import Path
except ImportError:
    print("Need Python 3.6 or higher, aborting")
    sys.exit(1)

if sys.version_info.major < 3:
    error("Need Python 3.6+, aborting")
    sys.exit(1)

elif sys.version_info.minor < 6:
    error("Needs Python 3.6+, aborting")
    sys.exit(1)


import tempfile
import argparse
import subprocess as sp
import datetime
import json
import fnmatch
import logging
import hashlib

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


# Determine default staging area
default_staging = "/tmp/{0}-lcdb-wf-staging".format(os.getenv('USER'))


def debug(s):
    logging.debug(GRAY + s + RESET)


def info(s):
    logging.info(GREEN + s + RESET)


def warning(s):
    logging.warning(YELLOW + s + RESET)


def error(s):
    logging.error(RED + s + RESET)



usage = f"""
This script assists in the deployment of relevant code from the lcdb-wf
repository to a new deployment directory for running an analysis.

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

"""


# Notes on include/exclude patterns for rsync:
#
# - Excluding a directory excludes everything below it
# - Including a directory does not automatically include everything below it.
# - Use "dir/***" to include everything
# - Patterns with no / applies to basename
# - Patterns ending with / implies directories only
# - Patterns starting with / implies that the root is the source dir provided to rsync
# - * is anything but /
# - ** is any part of path, including /
always = {
    "include": [
        "wrappers/wrappers",
        "include",
        "lib",
        "env.yaml",
        "env-r.yaml",
        ".gitignore",
    ],
    "exclude": [
        "sra_sampletable.tsv",
        "/.buildkite/*",
        "/ci/*",
        "/.circleci/*",
        "/config/*",
        "/deploy.py",
        "/docs/***",
        "/include/AnnotationHubCache",
        "/lib/postprocess/__pycache__",
        "/lib/__pycache__",
        "*/.pytest_cache/*",
        "/README.md",
        "/test/***",
        "/.travis.yml",
        "/workflows/*/results",
        "/workflows/*/data",
        "/workflows/figures/*",
        "/workflows/*/references_data",
        "/workflows/*/references_dir",
        "/workflows/*/reports",
        "/workflows/*/references/config",
        "/workflows/rnaseq/downstream/final_clusters",
        "/workflows/rnaseq/downstream/*html",
        "/workflows/rnaseq/downstream/*log",
        "/workflows/rnaseq/downstream/rnaseq_cache",
        "/workflows/rnaseq/downstream/rnaseq_files",
        "/workflows/rnaseq/downstream/*.tsv*",
        "/workflows/*/run*_test.sh",
        "/workflows/*/run_downstream_test.sh",
        "/workflows/*/Snakefile.test",
        "/workflows/*/.snakemake",
        "/wrappers/demo/*",
        "/wrappers/test/*",
        "/wrappers/test_toy.py",
    ],
}

flavors = {
    "chipseq": ["workflows/chipseq/*", "workflows/references/*"],
    "rnaseq": ["workflows/rnaseq/*", "workflows/references/*"],
    "colocalization": ["workflows/colocalization/*"],
    "full": ["workflows/*"],
}


def filter_out_excluded(filenames, patterns_to_exclude):
    """
    Only return the subset of `filenames` that do not match any
    `patterns_to_exclude`.
    """
    keep = []
    for filename in filenames:
        if not any(
            fnmatch.fnmatch(filename, pattern) for pattern in patterns_to_exclude
        ):
            keep.append(filename)
    return keep


def filter_out_other_workflows(filenames, patterns_to_include, prefilter="workflows/*"):
    """
    Return subset of `filenames` that:

        - don't match `prefilter`
        - match prefilter AND match any of `patterns_to_include`
    """
    keep = []
    for filename in filenames:
        # If it doesn't match prefilter, we don't want to do anything about it,
        # so keep it and move on.
        if not fnmatch.fnmatch(filename, prefilter):
            keep.append(filename)
            continue
        # Otherwise, only keep it if it's in patterns_to_include.
        if any(fnmatch.fnmatch(filename, pattern) for pattern in patterns_to_include):
            keep.append(filename)
    return keep


def write_file_list(source):

    # We start with anything under version control, and then progressively filter
    # out other stuff.
    under_version_control = sorted(
        sp.check_output(
            ["git", "ls-tree", "-r", "HEAD", "--name-only"],
            universal_newlines=True,
            cwd=source,
        ).splitlines(False),
    )
    keep = filter_out_excluded(under_version_control, always["exclude"])
    keep = filter_out_other_workflows(keep, flavors[flavor])

    exclude = tempfile.NamedTemporaryFile(delete=False).name
    with open(exclude, "w") as fout:
        fout.write("\n".join(always["exclude"]))

    include = tempfile.NamedTemporaryFile(delete=False).name
    with open(include, "w") as fout:
        fout.write("\n\n")
        fout.write("\n".join(keep) + "\n")

    debug("List of files excluded: {exclude}".format(**locals()))
    debug("List of files included: {include}".format(**locals()))

    return include, exclude


def clone_repo(dest, branch="master"):

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

    if this_md5 != that_md5:
        full_here = Path(__file__).resolve()
        full_there = Path(dest) / "deploy.py"
        error(
            "Files {full_here} and {full_there} do not match! ".format(**locals()) +
            "The deploy script you are running appears to be out of date. "
            "Please get an updated copy from https://github.com/lcdb/lcdb-wf, perhaps "
            "with 'wget https://raw.githubusercontent.com/lcdb/lcdb-wf/master/deploy.py'"
        )
        sys.exit(1)


def rsync(include, exclude, source, dest, rsync_args):
    rsync = [
        "rsync",
        "--relative",
        rsync_args,
        "--files-from={}".format(include),
        "--exclude-from={}".format(exclude),
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


def build_envs(dest, conda_frontend="mamba"):
    """
    Build conda environments.

    Parameters
    ----------

    dest : str
        Destination path. This is the project directory (likely as specified on
        the command line with --dest) in which the env and env-r yaml files
        should already exist. Envs will be created in here.

    conda_frontend : 'mamba' | 'conda'
        Which front-end to use (terminology borrowed from Snakemake)
    """
    mapping = [
        ("./env", "env.yml"),
        ("./env-r", "env-r.yml"),
    ]
    for env, yml in mapping:
        info("Building environment " + os.path.join(dest, env))

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
        help="""Options are {0}. Default is full.""".format(list(flavors.keys())),
    )
    ap.add_argument(
        "--dest", help="""Destination directory in which to copy files""", required=True
    )

    ap.add_argument(
        "--clone",
        help=f"""Make a new clone to a staging area (at the location specified
        by --staging which defaults to {default_staging}) and deploy from
        there. Useful if using this script as a standalone tool. You can also
        use --branch to configure which branch to deploy from that clone."""
    )

    ap.add_argument(
        "--staging",
        default=default_staging,
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
        default="mamba",
    )
    ap.add_argument(
        "--rsync-args",
        help="Options for rsync when deploying to a new directory. Default is %(default)s.",
        default="-rlt"
    )

    args = ap.parse_args()
    dest = args.dest
    flavor = args.flavor

    if args.staging:
        source = args.staging
        clone_repo(args.staging, args.branch)
    else:
        source = Path(__file__).parent.resolve()

    include, exclude = write_file_list(source)
    rsync(include, exclude, source, dest, args.rsync_args)
    deployment_json(source, dest)

    if args.build_envs:
        build_envs(dest, conda_frontend=args.conda_frontend)

    warning("Deployment complete in {args.dest}".format(**locals()))
