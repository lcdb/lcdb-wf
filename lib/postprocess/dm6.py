from snakemake.shell import shell

def fasta_postprocess(origfn, newfn):
    shell(
          "gunzip -fc {origfn} "
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

def fb_annotation_postprocess(origfn, newfn):
        shell(
            "gunzip -c {origfn} "
            "| grep -v '## '"
            r"| grep -v '\\'"
            "| tail -n +2"
            "| sed 's/#//g'"
            r"| sed 's/(s)//g'"
            "| gzip -c > {newfn} "
            "&& rm {origfn}"
        )

def fb_synonym_postprocess(origfn, newfn):
        shell(
            "gunzip -c {origfn} "
            "| awk '{{if ($1 ~ /^FBgn/ || $1 ~ /^##primary/){{print $0}}}}'"
            "| sed 's/##//g'"
            "| gzip -c > {newfn} "
            "&& rm {origfn}"
        )
