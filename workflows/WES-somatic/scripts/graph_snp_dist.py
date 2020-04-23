import pandas as pd
import matplotlib.pyplot as plt

def graph_snp_dist(gtf_vcf_bed, somaticonly_gtf_vcf_bed, sample, outfilename, somaticonly_outfilename):
    """
    This function will make a histogram of the distribution of variants across
    all genes. It takes two vcfs as inputs, one full one and one that has only
    the calls marked 'SOMATIC' by the variant callers. It also takes a genes
    only bed file that has been converted to a gtf
    """
    xlabel = 'Number of Snps in Gene'
    ylabel = 'Genes'
    df = pd.read_csv(gtf_vcf_bed, sep='\t')
    plt.hist(df.iloc[:,3], 50)
    plt.title(sample + ' Distribution of Snps across Genes')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(outfilename)
    plt.close()
    df = pd.read_csv(somaticonly_gtf_vcf_bed, sep='\t')
    plt.hist(df.iloc[:,3], 50)
    plt.title(sample + ' Distribution of Somatic Snps across Genes')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(somaticonly_outfilename)
    plt.close()


graph_snp_dist(snakemake.input[0], snakemake.input[1], snakemake.wildcards.sample, snakemake.output[0], snakemake.output[1])
