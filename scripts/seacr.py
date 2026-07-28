import os
import glob
from snakemake import shell

log = snakemake.log_fmt_shell()
logfile = None
extra = snakemake.params.get('extra', '')

outdir, basebed = os.path.split(snakemake.output.bed)
label = snakemake.params.block['label']
extra = snakemake.params.block.get('extra', '')

def bam_to_pe_bedgraph(bam):
    fragments_bedgraph = os.path.join(os.path.dirname(snakemake.output.bed), os.path.basename(bam) + ".fragments.bedgraph")
    shell(
        "bamCoverage "
        "--bam {bam} "
        "-o {fragments_bedgraph} "
        "-p {snakemake.threads} "
        "--outFileFormat bedgraph "
        "--skipNonCoveredRegions "  # SEACR expects zeros to be removed from the bedgraph

        # the following match the params used in the bigwig rule
        "--ignoreDuplicates "
        "--minMappingQuality 20 "
        "--binSize 1 "
        "--extendReads 300 "
    )
    return fragments_bedgraph

ip_bedgraph = bam_to_pe_bedgraph(snakemake.input.ip[0])
control_bedgraph = bam_to_pe_bedgraph(snakemake.input.control[0])


cmds = "SEACR_1.3.sh {ip_bedgraph} {control_bedgraph} norm relaxed $(dirname {snakemake.output})/peaks"
shell(cmds + ' {log}')
output = os.path.join(os.path.dirname(snakemake.output[0]), "peaks.relaxed.bed")
shell('''awk '{{OFS="\t"; print $1,$2,$3,".",$4,"."}}' {output} > {snakemake.output.bed}''')

# TODO: implement the bedgraph creation as a separate rule that runs if any seacr are configured.
# TODO: seacr creates tmp files in the working directory, which clutter the
# main workflow dir. They do appear to be cleaned up afterwards., but it would
# be better to set the working directory to the peak-calling dir so that the
# clutter only happens there.
shell('rm {ip_bedgraph} {control_bedgraph}')
