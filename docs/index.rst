lcdb-wf
=======

`lcdb-wf` is a collection of Snakemake workflows for common high-throughput
sequencing analysis, along with associated infrastructure.

We recognize that every experiment has its own unique idiosyncracies. Rather
than provide a one-size-fits-all solution, we aim to provide a reasonable
starting point that users can modify for their own use.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   getting-started
   guide
   workflows
   references
   wrappers
   tests


Guiding principles
------------------

- Expect users to know Snakemake

- Expect users will edit and customize the workflows accordingly

- Rely in wrappers where possible

- Test wrappers individually

- Test workflows to ensure end-to-end functionality

- Add complexity to wrappers rather than snakefile

- Wrappers should expose relatively simple APIs.

- Run references workflow once per site; other workflows either point directly
  to the created files or `include:` the references workflow to trigger updates

- Rely on tab-delimited sampletables as much as possible

- It's easier to delete than to create: each workflow will have "the works" and
  can be trimmed down according to the particular experiment's needs.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
