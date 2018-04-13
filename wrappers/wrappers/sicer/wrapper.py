import tempfile
import os
import glob
from snakemake import shell

log = snakemake.log_fmt_shell()
logfile = None
redundancy_threshold = snakemake.params.block.get('redundancy_threshold', snakemake.params.get('redundancy_threshold'))
window_size = snakemake.params.block.get('window_size', snakemake.params.get('redundancy_threshold'))
fragment_size = snakemake.params.block.get('fragment_size', snakemake.params.get('fragment_size'))
effective_genome_fraction = snakemake.params.block.get('effective_genome_fraction', snakemake.params.block.get('reference_effective_genome_fraction'))
gap_size = snakemake.params.block.get('gap_size', snakemake.params.get('gap_size'))
fdr = snakemake.params.block.get('fdr', snakemake.params.get('fdr'))
genome_build = snakemake.params.block.get('genome_build', snakemake.params.block.get('reference_genome_build'))

if redundancy_threshold is None:
    raise ValueError("SICER requires the specification of a 'redundancy_threshold'")
if window_size is None:
    raise ValueError("SICER requires the specification of a 'window_size'")
if fragment_size is None:
    raise ValueError("SICER requires the specification of a 'fragment_size'")
if effective_genome_fraction is None:
    raise ValueError("SICER requires the specification of an 'effective_genome_fraction'")
if gap_size is None:
    raise ValueError("SICER requires the specification of a 'gap_size'")
if fdr is None:
    raise ValueError("SICER requires the specification of an 'fdr'")
if genome_build is None:
    raise ValueError("SICER requires the specification of a recognized genome build")

outdir, basebed = os.path.split(snakemake.output.bed)
label = snakemake.params.block['label']

tmpdir = tempfile.mkdtemp()
cwd = os.getcwd()

cmds = (
    'bamToBed -i {snakemake.input.ip} > {tmpdir}/ip.bed ; '
    'bamToBed -i {snakemake.input.control} > {tmpdir}/in.bed '
)

shell("which python")
shell("whereis SICER")

shell(cmds)

#sicerexec = glob.glob(os.path.join(shell(whereis SICER | sed 's/\/bin\/SICER.sh/\/share/g' | awk '{print $2}'`, 'sicer*/SICER.sh'))
#print("the composed sicer executable is at")
#print(sicerexec)

cmds = (
    """cd {tmpdir} && """
    """whereis SICER | sed "s/\\/bin\\/SICER.sh/\\/share\\/sicer*\\/SICER.sh/g" | awk '{{printf("%s ", $2)}}' > run_command.bash && """
    """echo "{tmpdir} ip.bed in.bed {tmpdir} {genome_build} {redundancy_threshold} {window_size} {fragment_size} {effective_genome_fraction} {gap_size} {fdr}" >> run_command.bash && """
    """bash run_command.bash """
#    ' > tmp.sicer.output 2> tmp.sicer.error '
    """&& cd {cwd}"""
)

shell(cmds)

resultsfile = glob.glob(os.path.join(tmpdir, '*-islands-summary-FDR*'))

if len(resultsfile) == 1:
    hit = resultsfile[0]
    basehit = os.path.basename(resultsfile[0])
elif len(resultsfile) > 1:
    raise ValueError("Multiple islands-summary-FDR files found in temporary working directory: " + str(os.listdir(tmpdir)))
else:
    raise ValueError("No islands-summary-FDR file found: " + str(os.listdir(tmpdir)))

# Fix the output file so that it conforms to UCSC guidelines
#shell("mv {tmpdir}/tmp.sicer.output {snakemake.output.bed}.sicer.output")
#shell("mv {tmpdir}/tmp.sicer.error {snakemake.output.bed}.sicer.error")

shell(
    "export LC_COLLATE=C; "
    """awk -F"\\t" -v lab={label} '{{printf("%s\\t%d\\t%d\\t%s_peak_%d\\t%d\\t.\\t%g\\t%g\\t%g\\n", $1, $2, $3-1, lab, NR, -10*log($6)/log(10), $7, -log($6)/log(10), -log($8)/log(10))}}' """
    "{hit} > {snakemake.output.bed}.tmp "
    "&& bedSort {snakemake.output.bed}.tmp {snakemake.output.bed}"
    "&& rm {snakemake.output.bed}.tmp && rm -Rf {tmpdir}"
)
