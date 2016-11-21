A streamlined version of LCDB workflows.

Guiding principles:

- most of the work is in wrappers. The `lcdb-wrapper-tests` repo is included as a submodule here.

- in the previous `lcdb-workflows`, the idea was to have everything in a single
config file. While great for end-users, it was too much overhead to write new
functionality in the workflows and then add all the other infrastructure in
order to expose that new functionality to the config file.

- references workflow should ideally be run once per site; other workflows can
either point directly to the created files or can `include:` the workflow to
trigger updates

- make heavy use of sampletables

- any generally-useful helper functions go in `lcdblib`. We may want to include
that as submodule, too.

- It is expected that a particular workflow will get substantially edited
before actual use. Operating under the assumption that it's easier to delete
than to create, each workflow will have "the works" and can be trimmed down
according to the particular experiment's needs.

- Workflows should have a `patterns` section at the top that lays out, in one
place, the output file patterns. Rules can use these to define output patterns.
Using the `fill_patterns` function, these patterns are "rendered" into
a `targets` dictionary. The contents of the `targets` directory are then used
for the `all` rule. This provides, all in one section of the Snakefile, a fair
amount of customization options. Patterns can be commented out, or targets can
be excluded from the `all` rule to fine-tune which rules will be run.

- A side effect of the `patterns` and `targets` design is that aggregation
rules become much easier to write. If all FastQC runs are under a `fastqc` key
in `targets`, they can all be used as input to a rule using
`flatten(targets['fastqc'])`. This is in contrast to, say, writing input
functions.
