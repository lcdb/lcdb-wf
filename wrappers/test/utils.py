"""
Stripped-down version of Snakemake's test framework.
"""

import sys
import os
from textwrap import dedent
import subprocess as sp
import tempfile
import hashlib
import urllib
import shutil
import shlex
import inspect

import pytest
from snakemake import snakemake
from snakemake.shell import shell
from snakemake.utils import makedirs


SCRIPTPATH = shutil.which('snakemake')

# test data url
URL = 'https://github.com/lcdb/lcdb-test-data/blob/add-chipseq/data/{}?raw=true'


def tmpdir_for_func(factory):
    caller = inspect.stack()[1][3]
    return str(factory.mktemp(caller))


def _download_file(fn, d):
    """
    Intended to be called from a pytest.fixture function.

    `fn` is a path to a file that is used to fill in `URL`. `d` is a tempdir
    likely created by the calling function to which the file will be
    downloaded.

    The path to the downloaded file is returned.
    """
    url = URL.format(fn)
    dest = os.path.join(d, fn)
    makedirs(os.path.dirname(dest))
    basename = os.path.basename(fn)
    shell('wget -q -O- {url} > {dest}')
    return dest


def dpath(path):
    "path relative to this file"
    return os.path.realpath(os.path.join(os.path.dirname(__file__), path))


def md5sum(filename):
    data = open(filename, 'rb').read()
    return hashlib.md5(data).hexdigest()


def run(path, snakefile, check=None, input_data_func=None, tmpdir=None, use_conda=False, **params):
    """
    Parameters
    ----------

    path : str
        Path to a wrapper directory.

    snakefile : str
        Contents of a snakefile. `dedent()` will be run on it.

    check : callable or None
        After running the snakefile on the input data, this function will be
        called while inside the directory. This function is where the actual
        tests (assertions etc) should be performed.

        If None, the snakefile will be run but no tests will be performed on
        the output.

    input_data_func : None | callable
        If not None, then this callable object will be called with
        a single argument corresponding to the temp directory. It will be
        called after the wrapper and test-case contents have been copied to the
        temp dir, but before the test is run. It is expected to create any data
        required in whatever directory structure is required.

    tmpdir : None or path

    """
    # store any tempdirs here for later deletion
    to_clean_up = []


    if tmpdir is None:
        tmpdir = tempfile.mkdtemp(prefix='.test', dir=os.path.abspath('.'))
    else:
        tmpdir = str(tmpdir)
    try:
        # copy over the wrapper
        wrapper_dir = os.path.join(tmpdir, 'wrapper')
        os.makedirs(wrapper_dir)
        cmds = (
            'find {} -maxdepth 1 -type f -print0 | xargs -0 cp -t {}'
            .format(shlex.quote(path), shlex.quote(wrapper_dir))
        )
        sp.call(cmds, shell=True)

        # write the snakefile, filling in the "wrapper" placeholder
        with open(os.path.join(tmpdir, 'Snakefile'), 'w') as fout:
            fout.write('shell.executable("/bin/bash")\n')
            fout.write(dedent(snakefile))

        # Create the input data
        input_data_func(tmpdir)

        success = snakemake(os.path.join(tmpdir, 'Snakefile'), workdir=tmpdir, stats='stats.txt',
                            snakemakepath=SCRIPTPATH, config={}, use_conda=use_conda, **params)
        assert success, 'expected successful execution'

        # Change to the tmpdir and run the test function
        if check is not None:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            check()
            os.chdir(cwd)

    finally:
        for t in to_clean_up:
            shutil.rmtree(t)
        #shutil.rmtree(tmpdir)


def symlink_in_tempdir(mapping):
    """
    Returns a function that can be used for the `input_data_func` to utils.run.

    `mapping` is a dict where keys are 'target' and values are 'linkname'.

    It will symlink the data downloaded by the fixture into the temp dir
    created for the test case.
    """
    def _wrapped(tmpdir):
        for k, v in mapping.items():
            _linkname = os.path.join(tmpdir, v)
            _target = k
            _linkdir = os.path.dirname(_linkname)
            shell('mkdir -p {_linkdir} && ln -s {_target} {_linkname}')
    return _wrapped


def rm(path):
    shutil.rmtree(path)
