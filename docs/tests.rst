Tests
=====

These workflows and wrappers undergo continuous integration testing on
Travis-CI. Every push to GitHub results in the test suite being run on example
data.

There are some limitations on the free tier of travis-ci that we work around.
The main one is the 50-minute build time limit. We split the tests into roughly
independent jobs, and run them in parallel. There is a documentation-building
job, the py.test suite, and running the workflows on example data.

Another limitation is the practical size of the test data. It would be nice to
run tests on real-world data sets, but this is impractical to do in
a continuous integration context, and would waste the community resources on
travis-ci. Instead, we run the workflows on small example data sets built from
published data. These data sets, along with the Snakefile that generates them,
can be found in the `lcdb-test-data <https://github.com/lcdb/lcdb-test-data>`_
repo.

The workflows use test config and sampletable files in the ``config``
directory. The wrappers have their own test suite; see :ref:`wrappers` for
details.
