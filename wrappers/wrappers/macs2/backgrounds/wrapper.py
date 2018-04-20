import os
import os.path
from snakemake import shell
import tempfile
import subprocess
import fcntl

log = snakemake.log_fmt_shell()
logfile = None

lock_timeout_interval_seconds = 5
lock_timeout_max_intervals = 50

def conditionally_generate_files(_local,
                                 _d_estimate,
                                 _filename,
                                 _path,
                                 _out_bdg,
                                 _is_control):
    # check and see if the bdg file already exists
    sample_tag = "IP";
    if _is_control:
        sample_tag = "input"
    bdg_param_filename = _path + "/" + sample_tag + "_d" + str(_d_estimate)
    if _local is not None:
        bdg_param_filename = bdg_param_filename + "_local" + str(_local)
    bdg_param_filename = bdg_param_filename + ".bdg"
    # interpret (_local is None) as an extension by d
    local_est = _local
    if _local is None:
        local_est = _d_estimate

    # if the desired bdg file doesn't exist
    if not os.path.isfile(bdg_param_filename):
        # Prevent race conditions by getting a lock on the entire damn directory.
        # Unfortunately, due to making this exist outside of Snakemake, there is the
        # very real chance that a failed rule will cause some sort of consistency issue.
        # I think this solution works... mostly. There might still be a very small chance of collisions.
        lock_file = open(_path + "/.directory_lock", "a")
        for timeout in range(lock_timeout_max_intervals):
            try:
                # acquire directory lock
                fcntl.flock(lock_file, fcntl.LOCK_EX)
                # check one more time if it exists
                if os.path.isfile(bdg_param_filename):
                    fcntl.flock(lock_file, fcntl.LOCK_UN)
                    lock_file.close()
                    break
                # make the file
                local_intermediate_file = tempfile.NamedTemporaryFile().name
                # if it's a control, extend bidirectionally
                if _is_control:
                    shell("macs2 pileup -i {0} -B --extsize {1} -o {2} {3}".format(_filename, int(local_est/2), local_intermediate_file, log))
                    # if d was used, normalization is not required so save time
                    if _local is None:
                        shell("mv {0} {1}".format(local_intermediate_file, bdg_param_filename))
                    else:
                        shell("macs2 bdgopt -i {0} -m multiply -p {1} -o {2} {3}".format(local_intermediate_file,
                                                                                         _d_estimate/local_est,
                                                                                         bdg_param_filename,
                                                                                         log))
                        shell("rm {0}".format(local_intermediate_file))
                else:
                    shell("macs2 pileup -i {0} --extsize {1} -o {2} {3}".format(_filename, int(local_est), bdg_param_filename, log))
                # release directory lock
                fcntl.flock(lock_file, fcntl.LOCK_UN)
                lock_file.close()
                break
            except IOError:
                # in this case, the directory is in use already. Chill out for a reasonable amount of time.
                print("sleeping!")
                time.sleep(lock_timeout_interval_seconds)
            except Exception as e:
                # something unrelated failed. Clean up after yourself and then cry for help.
                fcntl.flock(lock_file, fcntl.LOCK_UN)
                lock_file.close()
                raise e
        # assuming everything's coded correctly, this catches when the mutex acquisition times out.
        if not lock_file.closed:
            fcntl.flock(lock_file, fcntl.LOCK_UN)
            lock_file.close()
            raise ValueError("macs2/background unable to generate {0} in reasonable time frame".format(bdg_param_filename))
    # symlink the output to the correct bdg file
    shell("ln -s {0} {1}".format(os.path.abspath(bdg_param_filename), _out_bdg))




d_estimate_file = snakemake.params.get("d_estimate")
d_estimate_default_value = snakemake.params.get("d_estimate_default_value", "300")
read_length_file = snakemake.params.get("read_length")
ip_bam = snakemake.input.get("ip_bam")
ip_bam_length = None
control_bam = snakemake.input.get("ctrl_bam")
control_bam_length = None

output_background_bdg_unscaled = snakemake.output.get("background_bdg_unscaled")
output_ip_bdg_unscaled = snakemake.output.get("ip_bdg_unscaled")
output_background_bdg_scaled = snakemake.output.get("background_bdg_scaled")
output_ip_bdg_scaled = snakemake.output.get("ip_bdg_scaled")

if d_estimate_file is None:
    raise ValueError("macs2/backgrounds requires params.d_estimate")
if read_length_file is None:
    raise ValueError("macs2/backgrounds requires params.read_length")
if ip_bam is None:
    raise ValueError("macs2/backgrounds requires input.ip_bam")
ip_bam_length = int(subprocess.check_output(["samtools", "view", "-c", ip_bam[0]]).decode('utf-8'))
if control_bam is None:
    raise ValueError("macs2/backgrounds requires input.ctrl_bam")
control_bam_length = int(subprocess.check_output(["samtools", "view", "-c", control_bam[0]]).decode('utf-8'))
if output_background_bdg_unscaled is None:
    raise ValueError("macs2/backgrounds requires output.background_bdg_unscaled")
if output_ip_bdg_unscaled is None:
    raise ValueError("macs2/backgrounds requires output.ip_bdg_unscaled")
if output_background_bdg_scaled is None:
    raise ValueError("macs2/backgrounds requires output.background_bdg_scaled")
if output_ip_bdg_scaled is None:
    raise ValueError("macs2/backgrounds requires output.ip_bdg_scaled")

# this wrapper requires control extensions for both slocal and llocal. These can be cached.
slocal = snakemake.params.block.get('slocal')
llocal = snakemake.params.block.get('llocal')
mfold_upper = snakemake.params.block.get('mfold_upper')
mfold_lower = snakemake.params.block.get('mfold_lower')

if slocal is None or \
   llocal is None:
    raise ValueError("macs2_extended runs require slocal and llocal parameters to be specified")
if mfold_lower is None or \
   mfold_upper is None:
    raise ValueError("macs2_extended runs require mfold lower and upper parameters to be specified")
slocal = slocal[0]
llocal = llocal[0]
mfold_lower = mfold_lower[0]
mfold_upper = mfold_upper[0]

effective_genome_count = snakemake.params.block.get('effective_genome_count',
                         snakemake.params.block.get('reference_effective_genome_count'))

genome_count_flag = ''
if effective_genome_count is not None:
    genome_count_flag = ' -g ' + effective_genome_count + ' '

# Multiple processes may try to write d estimates simultaneously.
# In this case, there's no control for if the d estimates are being made for the same set of parameters
# by multiple processes at the same time. In this case, the file may end up littered with repeat entries
# generated the very first time that the set of parameters is used. But that's not really so bad a thing,
# as the estimates will be identical and the first seen can be used in the future. Just make sure the file
# doesn't get corrupted.

# determine if this parameter set has already been estimated for this file

def try_to_find_d(_dfilename, _mfoldlower, _mfoldupper):
    if os.path.isfile(_dfilename):
        with open(_dfilename, "r") as f:
            for line in f:
                parsed = line.split()
                if len(parsed) == 3 and \
                   int(parsed[0]) == int(_mfoldlower) and \
                   int(parsed[1]) == int(_mfoldupper):
                    return True
    return False

value_available = try_to_find_d(d_estimate_file, mfold_lower, mfold_upper)

if not value_available:
    # Acquire the lock very, very early.
    lock_file = open(d_estimate_file + ".lock", "a")
    # give this process a reasonable number of attempts
    for timeout in range(lock_timeout_max_intervals):
        try:
            # attempt to acquire a lock.
            fcntl.flock(lock_file, fcntl.LOCK_EX)
            # this lock process could catch numerous processes the first time around.
            # you don't want to then have each and every one of the buggers compute their
            # own d estimates if you can avoid it.
            # so immediately check if d has appeared in the interim
            value_available = try_to_find_d(d_estimate_file, mfold_lower, mfold_upper)
            if value_available:
                fcntl.flock(lock_file, fcntl.LOCK_UN)
                lock_file.close()
                break
            
            cmds = (
                'macs2 predictd '
                '-i {ip_bam} '
                '{genome_count_flag} '
                '-m {mfold_lower} {mfold_upper} '
                '-f BAM '
            )
            shell(cmds + ' {log}')
            inlog = log.split(' ')[-2]
            d_estimate_temporary_file = tempfile.NamedTemporaryFile().name
            cmds = ('awk \'/predicted fragment length is/ {{print "{mfold_lower} {mfold_upper} "$(NF-1)}}\' '
                    '< {inlog} > {d_estimate_temporary_file}')
            shell(cmds)
            with open(d_estimate_temporary_file, "r") as f:
                lines = f.readlines()
                # try to get access to d estimate file
                with open(d_estimate_file, "a") as d_file:
                    # if there was less than one matching line from the awk extraction,
                    # it's because there was an MFOLD failure. Fall back to default extsize.
                    if len(lines) < 1:
                        d_file.write("{0} {1} {2}\n".format(mfold_lower, mfold_upper, d_estimate_default_value))
                    else:
                        d_file.write("{0} {1} {2}\n".format(mfold_lower, mfold_upper, lines[0].strip()))
            cmds = ('awk \'/tag size is determined as/ {{print $(NF-1)}}\' < {inlog} > {read_length_file}')
            shell(cmds)
            shell("rm {0}".format(d_estimate_temporary_file))

            # release the d estimate lock
            fcntl.flock(lock_file, fcntl.LOCK_UN)
            lock_file.close()
            break
        except IOError:
            # in this case, the lock file is in use already. Chill out for a reasonable amount of time.
            print("sleeping!!")
            time.sleep(lock_timeout_interval_seconds)
        except Exception as e:
            # something unrelated failed. Clean up after yourself and then cry for help.
            fcntl.flock(lock_file, fcntl.LOCK_UN)
            lock_file.close()
            raise e
    # assuming everything's coded correctly, this catches when the mutex acquisition times out.
    if not lock_file.closed:
        fcntl.flock(lock_file, fcntl.LOCK_UN)
        lock_file.close()
        raise ValueError("macs2/background unable to access d estimate file in reasonable time frame")


                

d_estimate = None
with open(d_estimate_file, "r") as f:
    for line in f:
        tokens = line.strip().split()
        if len(tokens) == 3 and \
           int(tokens[0]) == int(mfold_lower) and \
           int(tokens[1]) == int(mfold_upper):
            d_estimate = int(tokens[2])
            break
if d_estimate is None:
    raise ValueError("macs2/backgrounds requires preexisting d estimate for these parameters")


genome_background = float(control_bam_length) * float(d_estimate) / float(effective_genome_count)


slocal_llocal_merged = tempfile.NamedTemporaryFile().name
d_slocal_llocal_merged = tempfile.NamedTemporaryFile().name
slocal_background = tempfile.NamedTemporaryFile().name
llocal_background = tempfile.NamedTemporaryFile().name
d_background = tempfile.NamedTemporaryFile().name

macs2_ip_repo = os.path.split(ip_bam[0])[0] + "/macs2_bdgs/"
macs2_control_repo = os.path.split(control_bam[0])[0] + "/macs2_bdgs/"
shell("mkdir -p {0}".format(macs2_ip_repo))
shell("mkdir -p {0}".format(macs2_control_repo))

conditionally_generate_files(None,
                             d_estimate,
                             ip_bam,
                             macs2_ip_repo,
                             output_ip_bdg_unscaled,
                             False)
if not os.path.isfile(output_ip_bdg_unscaled):
    raise ValueError("unable to detect file output from ipbam extension")

conditionally_generate_files(None,
                             d_estimate,
                             control_bam,
                             macs2_control_repo,
                             d_background,
                             True)
if not os.path.isfile(d_background):
    raise ValueError("unable to detect file output from dbackground extension")

conditionally_generate_files(slocal,
                             d_estimate,
                             control_bam,
                             macs2_control_repo,
                             slocal_background,
                             True)
if not os.path.isfile(slocal_background):
    raise ValueError("unable to detect file output from slocal background extension")
conditionally_generate_files(llocal,
                             d_estimate,
                             control_bam,
                             macs2_control_repo,
                             llocal_background,
                             True)
if not os.path.isfile(llocal_background):
    raise ValueError("unable to detect file output from llocal background extension")
cmds = (
    'macs2 bdgcmp -m max '
    '-t {slocal_background} '
    '-c {llocal_background} '
    '-o {slocal_llocal_merged} && '
    'macs2 bdgcmp -m max '
    '-t {slocal_llocal_merged} '
    '-c {d_background} '
    '-o {d_slocal_llocal_merged} && '
    'macs2 bdgopt -m max '
    '-i {d_slocal_llocal_merged} '
    '-m max '
    '-p {genome_background} '
    '-o {output_background_bdg_unscaled} '
)

shell(cmds + '{log}')
shell("rm {0} {1} {2} {3} {4}".format(slocal_llocal_merged,
                                      d_slocal_llocal_merged,
                                      slocal_background,
                                      llocal_background,
                                      d_background))

if control_bam_length > ip_bam_length:
    mult_factor = ip_bam_length / control_bam_length
    shell('ln -s {0} {1}'.format(os.path.abspath(output_ip_bdg_unscaled),
                                 output_ip_bdg_scaled))
    cmds = (
        'macs2 bdgopt '
        '-i {output_background_bdg_unscaled} '
        '-m multiply '
        '-p {mult_factor} '
        '-o {output_background_bdg_scaled} ')
else:
    mult_factor = control_bam_length / ip_bam_length
    shell('ln -s {0} {1}'.format(os.path.abspath(output_background_bdg_unscaled),
                                 output_background_bdg_scaled))
    cmds = (
        'macs2 bdgopt '
        '-i {output_ip_bdg_unscaled} '
        '-m multiply '
        '-p {mult_factor} '
        '-o {output_ip_bdg_scaled} ')
shell(cmds + '{log}')
