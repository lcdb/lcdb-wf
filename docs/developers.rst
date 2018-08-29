Developers
==========

Adding a new aligner
--------------------

Modules
~~~~~~~

- in `lib/common.py`, there is a function `references_dict`. Within that is
  a `index_extensions` dictionary. You'll need to add the name of the aligner
  and the extension of the index it creates. If it creates multiple index
  files, just one should be sufficient. The filename will be automatically
  created and will be used as the expected output file which can then be
  accessed from the references dict as
  `references_dict[organism][tag][aligner]` for use in various rules that need
  the index as input (that is, any mapping rules).

Configuration
~~~~~~~~~~~~~

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
~~~~~~~

- For testing, create a copy of the config for any workflows it is used for,
  and change only the aligner.

- Modify `.circleci/config.yml` to include a new block in each of the
  variables, jobs, and workflows sections. Use the `rnaseq-star` blocks
  a guide for this. The idea is to only run up through the aligner step in
  a parallel task (to save on CI build time).
