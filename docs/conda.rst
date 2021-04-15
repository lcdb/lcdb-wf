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

You'll also probably want `mamba <https://github.com/mamba-org/mamba>`_. Mamba
is a drop-in replacement for conda that is faster and more robust. In fact, it
is now the default conda front-end for Snakemake. If you don't already have
mamba, you can install it into your base conda environment with:

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


If you use the ``--build-envs`` argument when deploying lcdb-wf to a project
directory (see :ref:`setup-proj`), two conda environments will be built in the
directories: ``env``, which has all of the non-R requirements, and ``env-r``
which has the R packages used in particular for downstream RNA-seq analysis.
These environments will use the fully-pinned environments in ``env.yml`` and
``env-r.yml``. If you've already deployed but didn't use the ``--build-envs``
argument, then then the equivalent command to run in the deployed directory is:

.. code-block:: bash

    mamba env create -p ./env --file env.yml
    mamba env create -p ./env-r --file env-r.yml

Troubleshooting environments
----------------------------

Sometimes there is a problem with creating an environment. For example, the
exact package specified in the env yaml might not be available for some reason
(this should not happen, but in practice sometimes it does in corner cases).

If this happens, you can try a couple things.

First, some terminology with how packages are specified in the environment
yamls. Here's an example for ``libpng`` version 1.6.37::

    libpng=1.6.37=hed695b0_2
    |____| |____| |________|
     |       |       |
    name     |       |
           version   |
                   build

The package name and version are pretty standard. The `build` string refers to
different built versions of the conda package, but for the same version (1.6.37
in this case) of the package. In this example, ``hed695b0`` is a hash of all
the pinned dependencies for this package at packaging time. The `conda-forge
pinning docs <https://conda-forge.org/docs/maintainer/pinning_deps.html>`_ give
more detail on what this pinning is about, but basically if that pinning
changes then this hash will change. The ``_2`` on the end of the build string
hash indicates that this is the third built package (build numbers start at
zero) for this version of ``libpng`` using the same pinning. In other words,
there also likely exists ``libpng=1.6.37=hed695b0_1`` and
``libpng=1.6.37=hed695b0_0``. At the time of this writing, there is also
``libpng-1.6.37-h21135ba_2`` (notice the different hash) which is the same
libpng version but uses different pinnings.

What does this mean for troubleshooting?

For any package that seems to be problematic, edit the respective environment
yaml (e.g., ``env.yml``) to remove the build string (so in the example above,
try changing to ``libpng=1.6.37``) and try building the environment again. If
that doesn't work, try removing the version as well (so just ``libpng``).

Alternatively for very problematic cases or cases where there are multiple
problematic packages, you can try creating an environment with the "loose"
pinning in ``include/requirements.txt`` which effectively does not require any
particular versions with the exception of a few corner cases. Keep in mind that
using that file may cause the environment to take a while to build as conda (or
mamba) solves the dependencies of all the specified packages.


Conda envs in lcdb-wf
---------------------

Given all of the software used across all of `lcdb-wf`, the environments can
take a lot of time to build because the solver needs to figure out the entire
dependency tree and come up with a solution that works to satisfy the entire
set of specified requirements.

We chose to split the conda environments in two: the **main** environment and the **R**
environment (see :ref:`conda-design-decisons`). These environments are
described by both "strict" and "loose" files. By default we use the "strict"
version, which pins all versions of all packages exactly. This is preferred
wherever possible. However we also provide a "loose" version that is not
specific about versions. The following table describes these files:

| strict version | loose version                  | used for                         |
+================+================================+==================================+
| ``env.yml``    | ``include/requirements.txt``   | Main Snakefiles                  |
| ``env-r.yaml`` | ``include/requirements-r.txt`` | Downstream RNA-seq analysis in R |

When deploying new instances, use the ``--build-envs`` argument which will use
the strict version. Or use the following commands in a deployed directory:

.. code-block:: bash

    mamba env create -p ./env --file env.yml
    mamba env create -p ./env-r --file env-r.yml

When getting ready to release a new lcdb-wf version, create a new environment
using the loose version to prepare the env and then when tests pass, export it
to yaml. That is:

.. code-block:: bash

    # use loose version when preparing a new version of lcdb-wf
    mamba create -p ./env --file include/requirements.txt
    mamba create -p ./env-r --file include/requirements-r.txt

    # then do testing....

    # when tests pass, export the envs
    conda env export -p ./env > env.yml
    conda env export -p ./env-r > env-r.yaml

    # commit, push, finalize release


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

.. note::

    Prior to v1.7, we used requirements.txt files with loose pinning. Moving to
    yaml files allows us the option of also installing pip packages if needed.
    It also allows us to specify channels directly in the yaml file for
    streamlined installation.

    Using strictly-pinned yaml files that are consistently tested will
    hopefully result in a more stable experience for users. For example, if you
    happen to create an environment around the time of a new R/Bioconductor
    release, the environment may not build correctly using a loose pinning.
    Other transient issues in the packaging ecosystem can similarly cause
    issues.
