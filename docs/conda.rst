.. _conda-envs:

Using conda and conda envs
==========================

Conda basics
------------

If you're not familiar with ``conda``, it installs particular versions of
software in an isolated location on your computer. When you "activate" the
environment, it places that location at the beginning of your ``$PATH``
variable, so that any executables there are found first. It does not affect
any existing installation of any software on your machine and does not need
root privileges.

If you don't already have conda installed and the Bioconda channel set up, see
the `Bioconda docs <https://bioconda.github.io>`_ for details.

**It is recommended that you create a separate environment directory for
each project**. That way you can update packages in each project
independently of any others, and yet the environment will always be close at
hand. This is an especially good practice in shared space as others can easily
find and activate the environment specific to the project.

Conda envs in lcdb-wf
---------------------

lcdb-wf has lot of requirements and it would be a lot of effort to install them
all. Furthermore, as you work on many projects over time, older projects may
need older versions, but you may want newer versions for more recent projects.

Conda environments handle all of this. Each project has its own isolated set of
software that is independent of other projects.

However, given all of the software used across all of lcdb-wf, the environments
can take a lot of time to build. Conda has to solve the entire dependency tree
and come up with a solution that works to satisfy the entire set of specified
requirements. That is, if we specify snakemake and STAR as requirements,
conda has to identify all of the dependencies of snakemake, all of the
dependencies of STAR, and inspect them to make sure they are compatible with
each other. If not, it has to incrementally try different versions of those
dependencies to find a solution. This process becomes exponentially more
complex as we add more requirementst. As a result, solving and building the
environments can take some time (usually under 15 mins).

The reason for having two different environments -- one for R, one for non-R --
is that by removing the entire sub-DAG of R packages from the main environment
we can dramatically reduce the creation time for an environment.

The **main** environment's requirements are stored in
``requirements-non-r.txt``. These are used for the primary workflows. The
**R-related** requirements are in ``requirements-r.txt`` and these are used for
downstream work, like ``workflows/rnaseq/downstream/rnaseq.Rmd`` and
``workflows/chipseq/downstream/diffbind.Rmd``.

Building the environments
-------------------------
If you use the ``--build-envs`` argument when deploying lcdb-wf to a project
directory (see :ref:`setup-proj`), conda environments will be built in the
directories ``env`` (with all non-R requirements) and ``env-r`` (R packages).

Otherwise, do the following in the top-level directory of the deployment:

.. code-block:: bash

    conda create -p ./env --file requirements-non-r.txt
    conda create -p ./env-r --file requirements-r.txt
