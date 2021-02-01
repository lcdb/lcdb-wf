For Developers
==============

.. _wrappers:

Wrappers
~~~~~~~~
One of the guiding principles of `lcdb-wf` is to put the complexity into
wrappers.  A `snakemake wrapper
<http://snakemake.readthedocs.io/en/latest/snakefiles/modularization.html#wrappers>`_
is a directory containing a script and a conda environment YAML.

The `wrappers/wrappers` directory contains a collection of wrappers that are
used in the workflows, or can be used in other workflows.

Writing wrappers
----------------
See the `demo wrapper
<https://github.com/lcdb/lcdb-wf/tree/master/wrappers/wrappers/demo>`_ for
a template to use when creating your own wrappers.

Running tests
-------------

While wrappers run in their own conda environment, we still have some
dependencies for the tests themselves. For example, one of the HISAT2 tests
verifies that the index has the correct chromosomes by running hisat2-inspect,
so hisat2 is a dependency.

With the bioconda channel set up, install the dependencies::

    conda create -n testenv --file requirements.txt

And run tests::

    source activate testenv
    pytest test/tests.py -v

Writing tests
-------------

Here, a *test* is a function that specifically checks output to ensure it is as expected.

A *fixture* is setup for a test. It can possibly itself be a test (if it checks
expected output) but doesn't have to be. An example of a fixture would be
a function that downloads a fastq file for use in other tests. Another example
of a fixture would be a function that runs cutadapt on that fastq, and provides
that trimmed fastq for other tests. An example of a test would be a function
that checks the cutadapt output to ensure it is as expected.

See ``test/test_demo.py`` for a heavily-annotated example, and
``test/raw_data_fixtures.py`` for useful fixtures.

Write new tests in ``tests/test_<wrappername>.py``. Put any fixtures at the top
of the file. Import any fixtures from other files as needed. There is
a ``raw_data_fixtures.py`` file that will be helpful.

A test consists of a *Snakefile* (typically written as a triple-quoted string
inside each test function), an *input data function* (typically created via
``utils.symlink_in_tempdir``), and a *test function* that performs arbitrary
tests on the output.

Importantly, for any test, you can go to
``/tmp/pytest-of-$USER/pytest-$N/$TEST_NAME`` directory and you'll see the
Snakefile, input data, and any created output files. This is useful for writing
the tests in the first place (e.g., what kinds of things to test for) or for
troubleshooting (e.g., go to the directory and run snakemake to see what
happens).

Some notes on how the fixtures work
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The py.test docs are not that great, so here's some extra info I've learned
that might be helpful for understanding how the tests work.

py.test fixtures are functions. They are used for things like setup and
teardown that needs to be done once, or for preparing data just once that can
be used across many tests. In this repo, the fixtures are primarily used for
downloading data from the testing data repo to a test-specific temporary
directory.

Fixtures are stored in a conftest.py file (see ``tests/conftest.py``). Any
discovered pytest tests in any module (functions starting with ``test_``) can
use any of these fixtures as an argument, and will get the return value of that
fixture at runtime.

pytest also has some built-in fixtures. Here we use the ``tmpdir`` fixture and
the ``tmpdir_factory`` fixture. The nice thing about these is that pytest
manages them on disk and only keeps the last handful of tempdirs. This makes it
straightforward to find the test data on disk and poke around when
troubleshooting.

Fixtures can use other fixtures, and what's really nice is that py.test handles
all of the dependencies across fixtures. For example, we can have a fixture
that downloads a FASTQ file and another fixture that runs FastQC on the FASTQ.
Then we can have a MultiQC test that uses the FastQC fixture as input. py.test
will figure that out and create and clean up tempdirs as needed.

I'm still working out the best way to do this. As one extreme, we could
maintain a single big Snakefile for each set of mutually compatible wrappers,
and returns a big dict of all the files created. This does a high-level
functional test of all of the wrappers, while the more specific wrapper tests
could verify the contents. The other extreme is to have a separate fixture for
every stage of the workflow.

We could also effectively reproduce the dependency DAG implied by a Snakefile
as py.test fixture dependencies. Not sure we need to go that route though.

Module documentation
~~~~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 2

   lib.common
   lib.chipseq
   lib.patterns_targets

Miscellanous tips
~~~~~~~~~~~~~~~~~

Adding a new aligner
--------------------

Modules
^^^^^^^

- in `lib/common.py`, there is a function `references_dict`. Within that is
  a `index_extensions` dictionary. You'll need to add the name of the aligner
  and the extension of the index it creates. If it creates multiple index
  files, just one should be sufficient. The filename will be automatically
  created and will be used as the expected output file which can then be
  accessed from the references dict as
  `references_dict[organism][tag][aligner]` for use in various rules that need
  the index as input (that is, any mapping rules).

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


