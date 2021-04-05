#!/usr/bin/env python
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
        "requirements-non-r.txt",
        "requirements-r.txt",
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
        "/workflows/rnaseq/downstream/final_clusters",
        "/workflows/rnaseq/downstream/*html",
        "/workflows/rnaseq/downstream/*log",
        "/workflows/rnaseq/downstream/rnaseq_cache",
        "/workflows/rnaseq/downstream/rnaseq_files",
        "/workflows/rnaseq/downstream/*.tsv*",
        "/workflows/*/run_test.sh",
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

ap = argparse.ArgumentParser(usage=usage)
ap.add_argument(
    "--flavor",
    default="full",
    help="""Options are {0}. Default is full.""".format(list(flavors.keys())),
)
ap.add_argument("--dest", help="""Destination directory in which to copy files""", required=True)
ap.add_argument(
    "--build-envs",
    action="store_true",
    help="""If specified, conda environments with all dependencies will be
    installed into directories called "env" and "env-r" within the directory
    provided for --dest.""",
)
ap.add_argument("--verbose", "-v", action="store_true", help="""Verbose mode""")

args = ap.parse_args()
dest = args.dest
flavor = args.flavor


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


# We start with anything under version control, and then progressively filter
# out other stuff.
under_version_control = sorted(
    sp.check_output(
        ["git", "ls-tree", "-r", "HEAD", "--name-only"], universal_newlines=True
    ).splitlines(False)
)
keep = filter_out_excluded(under_version_control, always["exclude"])
keep = filter_out_other_workflows(keep, flavors[flavor])

exclude = tempfile.NamedTemporaryFile(delete=False).name
if args.verbose:
    print("Exclude file: {}".format(exclude))
with open(exclude, "w") as fout:
    fout.write("\n".join(always["exclude"]))

include = tempfile.NamedTemporaryFile(delete=False).name
if args.verbose:
    print("Include file: {}".format(include))
with open(include, "w") as fout:
    fout.write("\n\n")
    fout.write("\n".join(keep) + "\n")

rsync = [
    "rsync",
    "--relative",
    "-ar",
    "--progress",
    "--files-from={}".format(include),
    "--exclude-from={}".format(exclude),
    HERE,
    dest,
]
if args.verbose:
    rsync.append("-vv")

sp.check_call(rsync)

# This next section builds the .lcdb-wf-deployment.json data.
#
# First, the last commit:
commit, message = (
    sp.check_output(["git", "log", "--oneline", "-1"], universal_newlines=True)
    .strip()
    .split(" ", 1)
)

# When we're deploying:
now = datetime.datetime.strftime(datetime.datetime.now(), "%Y%m%d%H%M")

# Where the remote was:
remotes = sp.check_output(["git", "remote", "-v"], universal_newlines=True)
remotes = [i.strip() for i in remotes.splitlines()]

# The branch we're deploying from:
branch = sp.check_output(["git", "branch"], universal_newlines=True)
branch = [i for i in branch.splitlines() if i.startswith("*")]
assert len(branch) == 1
branch = branch[0]
branch = branch.split("* ")[1]

d = {
    "git": {"commit": commit, "message": message, "remotes": remotes, "branch": branch},
    "timestamp": now,
}
log = os.path.join(dest, ".lcdb-wf-deployment.json")
with open(log, "w") as fout:
    fout.write(json.dumps(d) + "\n")
os.chmod(log, 0o440)


# If specified, build an environment in `dest/env`, using the correct channels.
if args.build_envs:
    sp.check_call(
        [
            "mamba",
            "env",
            "create",
            "-p",
            "./env",
            "--file",
            "env.yml",
            "-c",
            "conda-forge",
            "-c",
            "bioconda",
            "-c",
            "defaults",
        ],
        universal_newlines=True,
        cwd=dest,
    )
    sp.check_call(
        [
            "mamba",
            "env",
            "create",
            "-p",
            "./env-r",
            "--file",
            "env-r.yml",
            "-c",
            "conda-forge",
            "-c",
            "bioconda",
            "-c",
            "defaults",
        ],
        universal_newlines=True,
        cwd=dest,
    )
