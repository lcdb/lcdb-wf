import tempfile
import os
import glob
from snakemake import shell

logfile = None

# as SICER's interface is rather strict, this wrapper enforces named variables
# instead of 'extra' arbitrary string

def get_value(key, key2=None):
    """
    Get the value from params.block if it exists, otherwise from params.

    If key2 is not None, it's a different key to extract from the same params.block.

    Raises ValueError if nothing is configured.
    """
    if key2 is None:
        key2 = key
        val = snakemake.params.block.get(key, snakemake.params.get(key))
    else:
        val = snakemake.params.block.get(key, snakemake.params.block.get(key2))

    if val is None:
        raise ValueError(
            "SICER requires the specification of '{0}'".format(key))
    return val

redundancy_threshold = get_value('redundancy_threshold')
window_size = get_value('window_size')
fragment_size = get_value('fragment_size')
effective_genome_fraction = get_value('effective_genome_fraction', 'reference_effective_genome_fraction')
gap_size = get_value('gap_size')
fdr = get_value('fdr')
genome_build = get_value('genome_build', 'reference_genome_build')

outdir, basebed = os.path.split(snakemake.output.bed)
label = snakemake.params.block['label']

tmpdir = tempfile.mkdtemp()
cwd = os.getcwd()

# SICER expects bed input format, not bam as in other peak callers
shell(
    'bamToBed -i {snakemake.input.ip} > {tmpdir}/ip.bed ; '
    'bamToBed -i {snakemake.input.control} > {tmpdir}/in.bed '
)

# SICER emits a single hard-coded file that does not respect output directory.
# So move each run into its own temp directory to avoid collisions with
# other processes.
os.chdir(tmpdir)

shell(
    # there is a CI-specific bug, in which the python symlink is not correctly resolved to python2.7;
    # so as a really desperate hack, modify SICER's python calls to directly touch 2.7
    """sed 's/^python/$CONDA_PREFIX\/bin\/python2.7/' """
    """$CONDA_PREFIX/share/sicer*/SICER.sh > {tmpdir}/SICER.sh && chmod u+x {tmpdir}/SICER.sh """
)
shell(
    # run SICER
    """{tmpdir}/SICER.sh {tmpdir} ip.bed in.bed {tmpdir} """
    """{genome_build} {redundancy_threshold} {window_size} """
    """{fragment_size} {effective_genome_fraction} {gap_size} {fdr} > tmp.output 2>&1 """
)

# Move back once the run is complete.
os.chdir(cwd)

# one of the results files gets converted to the broadPeak format ala macs
resultsfile = glob.glob(os.path.join(tmpdir, '*-islands-summary-FDR*'))
if len(resultsfile) == 1:
    hit = resultsfile[0]
    basehit = os.path.basename(resultsfile[0])
elif len(resultsfile) > 1:
    raise ValueError(
        "Multiple islands-summary-FDR files found in {1}: {0}"
        .format(os.listdir(tmpdir), tmpdir)
    )
else:
    raise ValueError("No islands-summary-FDR file found in {1}: {0}".format(os.listdir(tmpdir), tmpdir))

# "summary graph for [the run] in bedGraph format"
summary_graph = glob.glob(os.path.join(tmpdir, '*-W{0}.graph*'.format(window_size)))
if len(summary_graph) == 1:
    summary_graph = summary_graph[0]
else:
    raise ValueError("SICER graph output file not found")

# the bedGraph file above, normalized by library size per million, in wig format
normalized_prefilter_wig = glob.glob(os.path.join(tmpdir, '*-W{0}-normalized.wig'.format(window_size)))
if len(normalized_prefilter_wig) == 1:
    normalized_prefilter_wig = normalized_prefilter_wig[0]
else:
    raise ValueError("SICER normalized prefilter wig file not found")

# "summary of all candidate islands with their statistical significance
candidate_islands = glob.glob(os.path.join(tmpdir, '*-W{0}-G{1}-islands-summary'.format(window_size, gap_size)))
if len(candidate_islands) == 1:
    candidate_islands = candidate_islands[0]
else:
    raise ValueError("SICER candidate islands file not found")

# "delineation of significant islands"
significant_islands = glob.glob(os.path.join(tmpdir, '*-W{0}-G{1}-FDR*-island.bed'.format(window_size, gap_size)))
if len(significant_islands) == 1:
    significant_islands = significant_islands[0]
else:
    raise ValueError("SICER significant islands file not found")

# "library of raw redundancy-removed reads on significant islands
redundancy_removed = glob.glob(os.path.join(tmpdir, '*-W{0}-G{1}-FDR*-islandfiltered.bed'.format(window_size, gap_size)))
if len(redundancy_removed) == 1:
    redundancy_removed = redundancy_removed[0]
else:
    raise ValueError("SICER redundancy removed library file not found")

# "wig file for the island-filtered redundancy-removed reads
normalized_postfilter_wig = glob.glob(os.path.join(tmpdir, '*-W{0}-G{1}-FDR*-islandfiltered-normalized.wig'.format(window_size, gap_size)))
if len(normalized_postfilter_wig) == 1:
    normalized_postfilter_wig = normalized_postfilter_wig[0]
else:
    raise ValueError("SICER normalized postfilter wig file not found")

shell(
    "export LC_COLLATE=C; "
    # format the output in broadPeak format
    # note that SICER can emit p-values of 0 and in that case this file will contain "inf" entries
    """awk -F"\\t" -v lab={label} """
    """'{{printf("%s\\t%d\\t%d\\t%s_peak_%d\\t%d\\t.\\t%g\\t%g\\t%g\\n", $1, """
    """$2, $3-1, lab, NR, -10*log($6)/log(10), $7, -log($6)/log(10), -log($8)/log(10))}}' """
    "{hit} > {snakemake.output.bed}.tmp && "
    # sort the bed file, just to be sure
    "bedSort {snakemake.output.bed}.tmp {snakemake.output.bed} && "
    # rename the assorted output files
    "mv {resultsfile} {snakemake.output.bed}-islands-summary-significant && "
    "mv {summary_graph} {snakemake.output.bed}.graph && "
    "wigToBigWig {normalized_prefilter_wig} {snakemake.input.chromsizes} {snakemake.output.bed}-normalized-prefilter.bigWig && "
    "wigToBigWig {normalized_postfilter_wig} {snakemake.input.chromsizes} {snakemake.output.bed}-normalized-postfilter.bigWig && "
    "mv {candidate_islands} {snakemake.output.bed}-islands-summary && "
    "mv {significant_islands} {snakemake.output.bed}-island.bed && "
    "mv {redundancy_removed} {snakemake.output.bed}-islandfiltered.bed && "
    "mv {tmpdir}/tmp.output {snakemake.output.bed}.log && "
    # clean up the temp directory
    "rm {snakemake.output.bed}.tmp && rm -Rf {tmpdir}"
)
