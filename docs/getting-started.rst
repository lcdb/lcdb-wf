.. _getting-started:

Initial setup
=============


Setup required once per system
------------------------------
The following needs to be performed on each system on which you will be running
the workflows.

1. Set up bioconda
~~~~~~~~~~~~~~~~~~

Follow the instructions for setting up `bioconda
<https://bioconda.github.io>`_.  This includes installing `conda` and setting
up the channels in the correct order.

This is required to be able to have all software automatically installed into
the working directory without needing admin rights on the machine.

2. Add the `lcdb` channel (one-time setup)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This enables the `lcdb` conda channel so that additional dependencies not
included in bioconda (primarily, the ``lcdblib`` package) can be installed:

.. code-block:: bash

    conda config --add channels lcdb


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


2. Create a new conda environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This creates a top-level environment with Snakemake and other requirements. It
should be `activated` any time you'll be working with these workflows.

If you're not familiar with ``conda``, it installs particular versions of
software in an isolated location on your computer. When you "activate" the
environment, it places that location at the beginning of your ``$PATH``
variable, so that any executables there are found first.

Here we're using the name "lcdb-wf" for the new environment, but you can use
anything. In fact, it's recommended that you choose a unique name for each
project. That way you can update packages in each project independently of any
others.

::

    conda create -n lcdb-wf -y python=3 --file requirements.txt

Then activate the environment::

    source activate lcdb-wf

Eventually when you're done, you can "deactivate", which removes the
environment location from your ``$PATH`` until the next time you activate it.
You might want to hold off on this for now::

    source deactivate

