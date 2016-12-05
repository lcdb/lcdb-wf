from snakemake.shell import shell

def fasta_postprocess(origfn, newfn):
    shell("""gunzip -c {origfn} | sed "s/>/>chr/g" | gzip -c > {newfn}  && rm {origfn}""")

def gtf_postprocess(origfn, newfn):
        shell(
            "gunzip -c {origfn} "
            """| awk -F "\\t" '{{OFS="\\t"; if ($8=="") $8="."; print "chr"$0}}' """
            "| gzip -c > {newfn} "
            "&& rm {origfn}")
