from snakemake.shell import shell

def fasta_postprocess(origfn, newfn):
    shell(
          "gunzip -c {origfn} "
          "| chrom_convert --from FlyBase --to UCSC --fileType FASTA -i - "
          "| gzip -c > {newfn} "
          "&& rm {origfn}")

def gtf_postprocess(origfn, newfn):
        shell(
            "gunzip -c {origfn} "
            """| awk -F "\\t" '{{OFS="\\t"; if ($8=="") $8="."; print $0}}' """
            """| awk -F "\\t" '{{OFS="\\t"; if ( ! ($7 != "-" && $9 ~ /FBgn0002781/)) print $0}}' """
            "| chrom_convert -i - --from FlyBase --to UCSC --fileType GTF"
            "| gzip -c > {newfn} "
            "&& rm {origfn}")
