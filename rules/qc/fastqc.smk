
# TODO: rm wrapper
rule fastqc:
    input:
        '{sample_dir}/{sample}/{sample}{suffix}'
    threads:
        6
    output:
        html='{sample_dir}/{sample}/fastqc/{sample}{suffix}_fastqc.html',
        zip='{sample_dir}/{sample}/fastqc/{sample}{suffix}_fastqc.zip',
    resources:
        mem_mb=gb(8),
        runtime=autobump(hours=2)
    script:
        utils.wrapper_for('fastqc/wrapper.py')
