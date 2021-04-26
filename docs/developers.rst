For Developers
==============

Creating and updating conda envs
--------------------------------

The ``env.yml`` and ``env-r.yml`` files contain fully-pinned versions of the
environments. This hopefully helps with stability and can dramatically speed up
the creation of environment. However these env definitions periodically need to
be updated.

To do so, create new environments using the unpinned versions in
``include/requirements.txt`` and ``include/requirements-r.txt``. This may take
substantially longer to create.

Then run the tests (:ref:`running-the-tests`) using those environments.

If all tests pass, then export the newly-created environments to the
``env.yml`` and ``env-r.yml`` files.

When you commit and push those files, the CI/CD system will detect that they
are different and will trigger a re-build of the cached environments and
proceed with the tests using those new environments.

Running the full complex datasets
---------------------------------

Prior to a release, the complex datasets should be run. These do a more
extensive job in testing the corner cases. This should be run on a cluster or
a machine with substantial resources. The configs can be found in
``include/test``. Here is how to run it using the WRAPPER_SLURM:

.. code-block:: bash

    sbatch ../../include/WRAPPER_SLURM \
      --configfile ../../test/test_configs/complex-dataset-rnaseq-config.yaml \
      --config sampletable=../../test/test_configs/complex-dataset-rnaseq-sampletable.tsv

Module documentation
--------------------

.. toctree::
   :maxdepth: 2

   lib.common
   lib.chipseq
   lib.patterns_targets


Adding a new aligner
--------------------

Modules
^^^^^^^

In `lib/common.py`, there is a function `references_dict`. Within that is
a `index_extensions` dictionary. You'll need to add the name of the aligner and
the extension of the index it creates. If it creates multiple index files, just
one should be sufficient. The filename will be automatically created and will
be used as the expected output file which can then be accessed from the
references dict as `references_dict[organism][tag][aligner]` for use in various
rules that need the index as input (that is, any mapping rules).

Configuration
^^^^^^^^^^^^^

- add the aligner to the `include/reference_configs/test.yaml` config file,
  "indexes:" section.

- write a rule in `workflows/references/Snakefile` to build the index. Use the
  other index-building rules there as a guide.

- Depending on which type of workflow the aligner is appropriate for, add
  a rule there. Enclose it in an "if:" clause to only run if the config file
  has specified that aligner.

- add the name to the list of supported aligners in `docs/config-yaml.rst`, in
  the "Aligner config" section.

- add appropriate memory/time requirements to `clusterconfig.yaml` for that
  aligner.

Testing
^^^^^^^

- For testing, create a copy of the config for any workflows it is used for,
  and change only the aligner.

- Modify `.circleci/config.yml` to include a new block in each of the
  variables, jobs, and workflows sections. Use the `rnaseq-star` blocks as
  a guide for this. The idea is to only run up through the aligner step in
  a parallel task (to save on CI build time).


.. _new-peak-caller:

Adding a new peak-caller
------------------------

First, write a wrapper for the peak-caller. You can use the ``macs2``, ``spp``,
and ``sicer`` wrappers as a guide. A wrapper should expect one or more sorted
and indexed BAM files as IP, one or more sorted and indexed BAM files as input.
The wrapper should create at least a sorted BED file of peaks, and can
optionally create other supplemental files as well.

Next, add the peak-caller to the top of ``lib/patterns_targets.py`` in the
``PEAK_CALLERS`` list.

Then write a rule for the peak-caller, again using ``macs2``, ``spp``, or
``sicer`` rules as a guide.

Last, add additional lines in
``workflows/chipseq/config/chipseq-patterns.yaml`` for the
``patterns_by_peaks`` key.

To test or use, add the new peak-caller to the
``workflows/chipseq/config/config.yaml`` file's ``peak_calling`` key.
