.. _conda-envs:

conda and conda envs in `lcdb-wf`
=================================

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

You'll also need `mamba <https://github.com/mamba-org/mamba>`_. Mamba is
a drop-in replacement for conda that is faster and more robust. In fact, it is
now the default conda front-end for Snakemake. If you don't already have mamba,
you can install it into your base conda environment with:

.. code-block:: bash

    conda install -n base -c conda-forge mamba


Building the environments
-------------------------

**It is recommended that you create a separate environment directory for
each project**, rather than a single environment for all projects. That way you
can update packages in each project independently of any others, and yet the
environment will always be close at hand. This is an especially good practice
in shared space as others can easily find and activate the environment specific
to the project.

.. note::

    We recommend using mamba rather than conda for the speed increase and
    ability to more correctly solve environments. See the `snakemake docs
    <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda>`_
    for more info.

We provide two different ways of specifying environments. One is a "strict"
pinning, which lists all installed versions of all dependencies (and
dependencies of dependencies...) to the exact version and build. This
configuration can be found in ``env.yml`` and ``env-r.yml``. The other
version is "relaxed", where the versions of tools are allowed to float to
whatever the latest version is at the time of environment creation.


If you use the ``--build-envs`` argument when deploying lcdb-wf to a project
directory (see :ref:`setup-proj`), conda environments will be built in the
directories ``env`` (with all non-R requirements) and ``env-r`` (R packages).
This will use the fully-pinned environments in ``env-pinned.yml`` and
``env-r-pinned.yml``. If you've already deployed without building the
environments, then the equivalent command is:

.. code-block:: bash

    mamba env create -p ./env --file env.yml
    mamba env create -p ./env-r --file env-r.yml

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

However, given all of the software used across all of `lcdb-wf`, the
environments can take a lot of time to build because the solver needs to figure
out the entire dependency tree and come up with a solution that works to
satisfy the entire set of specified requirements.

We chose to split the conda environments in two: the **main** environment and the **R**
environment (see :ref:`conda-design-decisons`). These environments are
described by both "strict" and "loose" files. By default we use the "strict"
version, which pins all versions of all packages exactly. This is preferred
wherever possible. However we also provide a "loose" version that is not
specific about versions. The following table describes these files:

| strict version | loose version          | used for                         |
+================+========================+==================================+
| ``env.yml``    | ``requirements.txt``   | Main Snakefiles                  |
| ``env-r.yaml`` | ``requirements-r.txt`` | Downstream RNA-seq analysis in R |

When working on a new version of `lcdb-wf`, we use the loose version to
generate a new environment using the latest versions of all software. When the
tests are confirmed to be passing, the definition of these environments is
exported to the strict version for use in production.


.. _conda-design-decisions:

Design decisions
----------------

We made the design decision to split the conda envs into two different
environments -- one for R, one for non-R. We round that by by removing the
entire sub-DAG of R packages from the main environment we can dramatically
reduce the creation time.

We also made the decision to use large top-level environments rather than
smaller environments created for each rule using the ``conda:`` directive. This
allows us to activate a single environment to give us access to all the tools
used. This streamlines troubleshooting because we don't have to dig through the
``.snakemake/conda`` directory to figure out which hash corresponds to which
file, but comes with the up-front cost of creating the environment initially.

Building the environments
-------------------------
If you use the ``--build-envs`` argument when deploying lcdb-wf to a project
directory (see :ref:`setup-proj`), conda environments will be built in the
directories ``env`` (with all non-R requirements) and ``env-r`` (R packages).

Otherwise, do the following in the top-level directory of the deployment:

.. code-block:: bash

    # if you don't already have mamba:
    conda install mamba -c conda-forge

    mamba env create -p ./env --file env.yml
    mamba env create -p ./env-r --file env-r.yml
