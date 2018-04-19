import os
import os.path
from snakemake import shell
import tempfile
import subprocess

log = snakemake.log_fmt_shell()
logfile = None
extra = snakemake.params.get('extra', '')

def conditionally_generate_files(_local,
                                 _mfold_lower,
                                 _mfold_upper,
                                 _d_estimate,
                                 _ctrl_filename,
                                 _path,
                                 _out_bdg):
    # find the appropriate estimate of d

    # check and see if the bdg file already exists
    bdg_param_filename = _path + "/mfold_lower" + str(_mfold_lower) + "_upper" + str(_mfold_upper) + "_d" + str(_d_estimate)
    if _local is not None:
        bdg_param_filename = bdg_param_filename + "_local" + str(_local)
    bdg_param_filename = bdg_param_filename + ".bdg"
    # interpret (_local is None) as an extension by d
    local_est = _local
    if _local is None:
        local_est = _d_estimate
    # if the file doesn't exist
    if not os.path.isfile(bdg_param_filename):
        # make the file
        tmp = tempfile.NamedTemporaryFile().name
        shell("macs2 pileup -i {0} -B --extsize {1} -o {2} {3}".format(_ctrl_filename, int(local_est/2), tmp, log))
        # if d was used, normalization is not required so save time
        if _local is None:
            shell("mv {0} {1}".format(tmp, bdg_param_filename))
        else:
            shell("macs2 bdgopt -i {0} -m multiply -p {1} -o {2} {3}".format(tmp,
                                                                             _d_estimate/local_est,
                                                                             bdg_param_filename,
                                                                             log))
            shell("rm {0}".format(tmp))
    # symlink the output to the correct bdg file
    shell("ln -s {0} {1}".format(os.path.abspath(bdg_param_filename), _out_bdg))


label = snakemake.params.block['label']
extra = snakemake.params.block.get('extra', '')

d_estimate_filename = snakemake.input.d_estimate
control_bam = snakemake.input.ctrl_bam
control_bam_length = None
output_bdg = snakemake.output.bdg

if d_estimate_filename is None:
    raise ValueError("macs2/combine_backgrounds requires input.d_estimate")
if control_bam is None:
    raise ValueError("macs2/combine_backgrounds requires input.ctrl_bam")
else:
    control_bam_length = subprocess.check_output(["samtools", "view", "-c", control_bam[0]]).decode('utf-8')
if output_bdg is None:
    raise ValueError("macs2/combine_backgrounds requires output.bdg")

# this wrapper requires control extensions for both slocal and llocal. These can be cached.
slocal = snakemake.params.block.get('slocal')[0]
llocal = snakemake.params.block.get('llocal')[0]
mfold_upper = snakemake.params.block.get('mfold_upper')[0]
mfold_lower = snakemake.params.block.get('mfold_lower')[0]

if slocal is None or \
   llocal is None:
   raise ValueError("macs2_extended runs require slocal and llocal parameters to be specified")


effective_genome_count = snakemake.params.block.get('effective_genome_count',
                         snakemake.params.block.get('reference_effective_genome_count', ''))



d_estimate = None
with open(d_estimate_filename) as f:
    for line in f:
        tokens = line.strip().split()
        if len(tokens) == 3 and \
           int(tokens[0]) == int(mfold_lower) and \
           int(tokens[1]) == int(mfold_upper):
            d_estimate = int(tokens[2])
            break
if d_estimate is None:
    raise ValueError("macs2/combine_backgrounds requires preexisting d estimate for these parameters")


genome_background = float(control_bam_length) * float(d_estimate) / float(effective_genome_count)


tmp1 = tempfile.NamedTemporaryFile().name
tmp2 = tempfile.NamedTemporaryFile().name
tmp3 = tempfile.NamedTemporaryFile().name
tmp4 = tempfile.NamedTemporaryFile().name
tmp5 = tempfile.NamedTemporaryFile().name


macs2_repo = os.path.split(control_bam[0])[0] + "/macs2_bdgs/"
shell("mkdir -p {0}".format(macs2_repo))
conditionally_generate_files(None,
                             mfold_lower,
                             mfold_upper,
                             d_estimate,
                             control_bam,
                             macs2_repo,
                             tmp5)

conditionally_generate_files(slocal,
                             mfold_lower,
                             mfold_upper,
                             d_estimate,
                             control_bam,
                             macs2_repo,
                             tmp3)

conditionally_generate_files(llocal,
                             mfold_lower,
                             mfold_upper,
                             d_estimate,
                             control_bam,
                             macs2_repo,
                             tmp4)

cmds = (
    'macs2 bdgcmp -m max '
    '-t {tmp3} '
    '-c {tmp4} '
    '-o {tmp1} && '
    'macs2 bdgcmp -m max '
    '-t {tmp1} '
    '-c {tmp5} '
    '-o {tmp2} && '
    'macs2 bdgopt -m max '
    '-i {tmp2} '
    '-m max '
    '-p {genome_background} '
    '-o {output_bdg} '
)
shell(cmds + '{log}')
shell("rm {0} {1} {2} {3} {4}".format(tmp1, tmp2, tmp3, tmp4, tmp5))
