import os
import tempfile
from snakemake.shell import shell
# Inspired by http://wresch.github.io/2014/01/31/merge-bigwig-files.html

# If memory was supplied, we'll use that for sorting.
if 'memory' in snakemake.params:
    mem_arg = '-S {snakemake.params.memory}'
else:
    mem_arg = ''

# Copied from lcdblib to avoid having it as a dependency just for this
# function.
def make_relative_symlink(target, linkname):
    linkdir = os.path.dirname(linkname)
    relative_target = os.path.relpath(target, start=linkdir)
    linkbase = os.path.basename(linkname)
    if not os.path.exists(linkdir):
        shell('mkdir -p {linkdir}')
        shell('cd {linkdir}; ln -sf {relative_target} {linkbase}')

if len(snakemake.input.bigwigs) == 1:
    make_relative_symlink(snakemake.input.bigwigs[0], snakemake.output[0])

else:

    # bigWigMerge outputs sum; we need to divide each by n.
    f = 1.0 / len(snakemake.input.bigwigs)

    tmp = tempfile.NamedTemporaryFile(delete=False).name
    tmpdir = tempfile.gettempdir()

    shell(
        'export LC_ALL=C; '
        'bigWigMerge {snakemake.input.bigwigs} stdout '
        """| awk 'BEGIN{{OFS="\t"}}{{$4={f}*$4; print}}' """
        '| sort {mem_arg} -T {tmpdir} -k1,1 -k2,2n > {tmp} '
        '&& bedGraphToBigWig {tmp} {snakemake.input.chromsizes} '
        '{snakemake.output} '
    )
