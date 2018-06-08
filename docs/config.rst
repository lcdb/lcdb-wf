
.. _config:


Configuration
=============

The majority of the work in setting up a new project is in the configuration --
which samples to run, where the data files are, which references are
needed, etc.

**The entry point for configuration** is in the ``config/config.yaml`` file found
in each workflow directory. See :ref:`config-yaml` for more.

.. toctree::
   :maxdepth: 2

   config-yaml

The **references section** of the config file configures the genomes,
transcriptomes, and annotations to be used. See :ref:`references-config` for more.

.. toctree::
   :maxdepth: 2

   references-config

The **sample table**, lists sample IDs, filenames, and other metadata. Its path
is specified in the config file. See :ref:`sampletable` for more.

.. toctree::
   :maxdepth: 2

   sampletable

A **patterns file** only needs to be edited
if you're doing custom work. It determines the patterns of files that will be
created by the workflow. Ssee :ref:`patterns-and-targets` for more.

.. toctree::
   :maxdepth: 2

   patterns-targets
