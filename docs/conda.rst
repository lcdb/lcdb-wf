.. _conda-envs:

Using conda and conda envs
==========================

Conda basics
------------

If you're not familiar with ``conda``, it is a way of keeping software isolated
on a computer in an "environment" (basically a directory with the executables
for all the software you want to use). When you "activate" the environment, it
places that location at the beginning of your ``$PATH`` variable, so that any
executables there are found first. It does not affect any existing installation
of any software on your machine and does not need root privileges.

If you don't already have conda installed and the Bioconda channel set up, see
the `Bioconda docs <https://bioconda.github.io>`_ for details.

**It is recommended that you create a separate environment directory for
each project**. That way you can update packages in each project
independently of any others, and yet the environment will always be close at
hand. This is an especially good practice in shared space as others can easily
find and activate the environment specific to the project.

Building the environments
-------------------------

.. note::

    We recommend using mamba rather than conda for the speed increase and
    ability to more correctly solve environments. See the `snakemake docs
    <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda>`_
    for more info.


We provide two different ways of specifying environments. One is a "strict"
pinning, which lists all installed versions of all dependencies (and
dependencies of dependencies...) to the exact version and build. The other
version is "relaxed", where the versions of tools are allowed to float to
whatever the latest version is at the time of environment creation.

If you use the ``--build-envs`` argument when deploying lcdb-wf to a project
directory (see :ref:`setup-proj`), conda environments will be built in the
directories ``env`` (with all non-R requirements) and ``env-r`` (R packages).
This will use the fully-pinned environments in ``env-pinned.yml`` and
``env-r-pinned.yml``. If you've already deployed without building the
environments, then the equivalent command is:

.. code-block:: bash

    mamba env create -p ./env --file env-pinned.yml
    mamba env create -p ./env-r --file env-r-pinned.yml

If you want to use the more relaxed environment specifications, then use the
non-pinned yaml files. We recommend using the pinned version though, since that
is what is used for the automated tests.

.. conda-block:: bash

    mamba env create -p ./env --file env.yml
    mamba env create -p ./env-r --file env-r-pinned.yml

.. note::

    Prior to v1.7, we used requirements.txt files with loose pinning. Moving to
    yaml files allows us the option of also installing pip packages if needed.
    It also allows us to specify channels directly in the yaml file for
    streamlined installation.

    In addition, encouraging the use of strictly-pinned yaml files that are
    consistently tested will hopefully result in a more stable experience for
    users. For example, if you happen to create an environment around the time
    of a new R/Bioconductor release, the environment may not build correctly.
    Other transient issues in the packaging ecosystem can similarly cause
    issues.


Conda envs in lcdb-wf
---------------------

lcdb-wf has lot of requirements and it would be a lot of effort to install them
all. Furthermore, as you work on many projects over time, older projects may
need older versions, but you may want newer versions for more recent projects.

Conda environments handle all of this. Each project has its own isolated set of
software that is independent of other projects.

It can take some of time to build an environment when there are many
dependencies like in the environments `lcdb-wf` uses. This is because conda has
to solve the entire dependency tree and come up with a solution that works to
satisfy the entire set of specified requirements. Say we specify snakemake and
STAR as requirements. Then conda has to identify all of the dependencies of
snakemake, all of the dependencies of STAR, and inspect them to make sure they
are compatible with each other. If not, it has to incrementally try different
versions of those dependencies to find a solution. This process becomes
exponentially more complex as we add more requirements. As a result, solving
and building the environments can take some time (usually under 15 mins).

The reason for having two different environments -- one for R, one for non-R --
is that by removing the entire sub-DAG of R packages from the main environment
we can dramatically reduce the creation time for an environment.

The **main** environment's requirements are stored in ``env.yml``, or the
strictly pinned version ``env-pinned.yml``.  These are used for the primary
workflows. The **R-related** requirements are in ``env-r.yml`` (or
``env-r-pinned.yml``) and these are used for downstream work, like
``workflows/rnaseq/downstream/rnaseq.Rmd`` and
``workflows/chipseq/downstream/diffbind.Rmd``.

