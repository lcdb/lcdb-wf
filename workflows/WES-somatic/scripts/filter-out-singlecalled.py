

def filter_out_singlecalled(vcf, samplename, outfilename):
    with open(outfilename, 'w') as fout:
        for line in open(vcf):
            if line[0] == '#':
                fout.write(line)
            else:
                isec = (
                    line.strip()
                    .split('\t')
                    [7].split('set=')
                    [1].split(';')[0]
                )
                if '-' in isec or 'Intersection' in isec:
                    fout.write(line)

filter_out_singlecalled(snakemake.input[0], snakemake.wildcards.sample, snakemake.output[0])

#if __name__ == '__main__':
#    import argparse
#    ap= argparse.ArgumentParser()
#    ap.add_argument('vcf'),
#    ap.add_argument('samplename')
#    args=ap.parse_args()
#    filter_out_singlecalled(
#        args.vcf,
#        args.samplename
#    )



                
