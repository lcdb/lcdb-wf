.. _getting-started:

Initial setup
=============

.. note::

    `lcdb-wf` is tested and heavily used on Linux.

    It is likely to work on macOS as long as all relevant conda packages are
    available for macOS -- though this is not tested.

    It will **not** work on Windows due to a general lack of support of Windows
    in bioinformatics tools.

Setup required once per system
------------------------------
We use `bioconda <https://bioconda.github.io>`_ to automatically install
software into the working directory without needing admin rights on the
machine.

If you have the Anaconda Python distribution, you already have conda.
Otherwise, install `Miniconda <https://conda.io/miniconda.html>`_.

**Optional:** 
- Run the following commands to set up your conda channels. This puts the
  `lcdb` channel as lowest priority (this channel has the `lcdblib` package),
  and then matches the channel order required by `bioconda`.

.. code-block:: bash

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge


Setup required once per project
-------------------------------

The following should be done each time you set up a project.

1. Clone the git repo
~~~~~~~~~~~~~~~~~~~~~

Clone the repository from github into a new directory of your choosing and
change to that directory.

.. code-block:: bash

    git clone https://github.com/lcdb/lcdb-wf.git my-project-dir
    cd my-project-dir


.. _create-env:

2. Create a new conda environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following command will create a top-level environment with Snakemake and
other requirements. It needs to be `activated` any time you'll be working with
these workflows.

If you're not familiar with ``conda``, it installs particular versions of
software in an isolated location on your computer. When you "activate" the
environment, it places that location at the beginning of your ``$PATH``
variable, so that any executables there are found first. It does not affect any
existing installation of any software on your machine and does not need root
privileges.

**It is recommended that you use a different environment name for each
project**. That way you can update packages in each project independently of
any others. Here we use the name "lcdb-wf" for the new environment, but you can
use anything. Note that here we specify the channels to use, which include
``bioconda`` which depends on ``conda-forge``, and ``lcdb`` which provides the
``lcdblib`` package used by these workflows.

::

    conda create -n lcdb-wf --file requirements.txt

Then activate the environment::

    source activate lcdb-wf

Eventually when you're done, you can "deactivate", which removes the
environment location from your ``$PATH`` until the next time you activate it.
You might want to hold off on this for now if you'll be running the tests::

    source deactivate

.. note::

   An alternative approach is to create an environment at a specific path, for
   example inside a project directory:

   .. code-block:: bash

       conda create -p ./env --file requirements.txt
       source activate ./env

Next steps
----------

You may want to run tests to make sure everything is set up (see
:ref:`running-the-tests`), or jump right in to learning about how to configure
the workflows for your particular experiment (see :ref:`config`).
