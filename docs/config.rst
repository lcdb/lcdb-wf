
.. _config:


Configuration
=============

The majority of the work in setting up a new project is in the configuration --
which samples to run, where the data files are, which references are
needed, etc.

The entry point for configuration is in the ``config/config.yaml`` file found
in each workflow directory (see :ref:`config-yaml`).

The config file has a references section (see :ref:`references-config`) to
configure the genomes, transcriptomes, and annotations to be used.

The config file also points to a sampletable (see :ref:`sampletable`) which
lists sample IDs, filenames, and other metadata.

A patterns file (see :ref:`patterns-and-targets`) determines the patterns of
files that will be created by the workflow (and doesn't normally need to be
changed).

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   config-yaml
   sampletable
   references-config
   patterns-targets

