import os
import tempfile


def tempdir_for_biowulf():
    """
    Get an appropriate tempdir.

    The NIH biowulf cluster allows nodes to have their own /lscratch dirs as
    local temp storage. However the particular dir depends on the slurm job ID,
    which is not known in advance. This function sets the shell.prefix to use
    such a tempdir if it's available; otherwise it leaves TMPDIR unchanged.
    This makes it suitable for running locally or on other clusters, however if
    you need different behavior for a different cluster, then a different
    function will need to be written.
    """
    tmpdir = tempfile.gettempdir()
    jobid = os.getenv('SLURM_JOBID')
    if jobid:
        tmpdir = os.path.join('/lscratch', jobid)
    return tmpdir
