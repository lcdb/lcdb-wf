Features common to all workflows
================================

`patterns` and `targets`
------------------------
Workflows have a `patterns` dict at the top that lays out, in one place, the
output file patterns. Rules can optionally use values from this dict to define
output patterns. Using the `lcdblib.snakemake.fill_patterns` function, these
patterns are "rendered" into a `targets` dictionary (this is basically
a recursive `expand()`). Selected contents of the `targets` directory are then
used for the `all` rule.

This provides, all in one section of the Snakefile, a fair amount of
customization options. Patterns can be commented out, or targets can be excluded
from the `all` rule to fine-tune which rules will be run.

A nice side-effect of the `patterns` and `targets` design is that aggregation
rules become much easier to write. If all FastQC runs are under a `fastqc` key
in `targets`, they can all be used as input to a rule using
`lcdblib.utils.flatten(targets['fastqc'])`. This is in contrast to, say, writing
input functions.
