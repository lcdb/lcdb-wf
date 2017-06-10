__author__ = "Ryan Dale"
__copyright__ = "Copyright 2016, Ryan Dale"
__email__ = "dalerr@niddk.nih.gov"
__license__ = "MIT"

from snakemake.shell import shell

extra = snakemake.params.get('extra', '')
log = snakemake.log_fmt_shell()
inputs = snakemake.input
outputs = snakemake.output

if isinstance(inputs, dict) and isinstance(outputs, dict):
    # Get inputs
    in_R1 = inputs.get('R1', None)
    in_R2 = inputs.get('R2', None)
    in_FASTQ = inputs.get('fastq', None)

    if (in_R1 is None) and (in_FASTQ is not None):
        in_R1 = in_FASTQ
    elif (in_R1 is None) and (in_FASTQ is None):
        raise KeyError('If providing a dictionary for input/output, you must uese either '
            '`R1` or `fastq` for the first read. If providing a second read you must use `R2`.')

    # Get outputs
    out_R1 = outputs.get('R1', None)
    out_R2 = outputs.get('R2', snakemake.params.get('R2', None))
    out_FASTQ = outputs.get('fastq', None)

    if (out_R1 is None) and (out_FASTQ is not None):
        out_R1 = out_FASTQ
    elif (out_R1 is None) and (out_FASTQ is None):
        raise KeyError('If providing a dictionary for input/output, you must uese either '
            '`R1` or `fastq` for the first read. If providing a second read you must use `R2`.')

elif isinstance(inputs, list) and isinstance(outputs, list):
    # Get inputs
    if len(inputs) == 1:
        in_R1 = inputs[0]
        in_R2 = None
    elif len(inputs) == 2:
        in_R1 = sorted(inputs)[0]
        in_R2 = sorted(inputs)[1]
    else:
        raise IndexError("If providing a list for input/output, they must have either 1 or 2 values.")

    # Get outputs
    if len(outputs) == 1:
        out_R1 = outputs[0]
        out_R2 = snakemake.params.get('R2', None)
    elif len(outputs) == 2:
        out_R1 = sorted(outputs)[0]
        out_R2 = sorted(outputs)[1]
    else:
        raise IndexError("If providing a list for input/output, they must have either 1 or 2 values.")

# Run paired end if both in_R2 and out_R2 are provided
if (in_R2 is not None) and (out_R2 is not None):
    shell(
        "atropos trim "
        "--threads {snakemake.threads} "
        "{extra} "
        "-pe1 {in_R1} "
        "-pe2 {in_R2} "
        "-o {out_R1} "
        "-p {out_R2} "
        "{log}"
    )
elif (in_R1 is not None) and (out_R1 is not None) and (in_R2 is None) and (out_R2 is None):
    shell(
        "atropos trim "
        "{extra} "
        "--threads {snakemake.threads} "
        "-se {in_R1} "
        "-o {out_R1} "
        "{log}"
    )
else:
    raise ValueError("Input and Output must match. If you give two value for "
        "input you must give two values for output.")
