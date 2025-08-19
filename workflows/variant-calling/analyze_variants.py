# Import libraries
from cyvcf2 import VCF, Writer
from collections import defaultdict
import pickle
import re
import pandas as pd
import numpy as np
import re
import sys
import argparse
import yaml
import os
import csv
import math

def main():
    parser = argparse.ArgumentParser(description="A script to compute metrics for variants")

    parser.add_argument("--vcf", "-v", default="results/annotated/filtered_ann.vcf.gz", help="Path to vcf file containing variant information")
    parser.add_argument("--output", "-o", default="results/analyzed", help="Path to the folder that the output files will be written to")
    parser.add_argument("--gnomad", "-g", action='store_true', default=True, help="Flag to set if vcf is annotated with gnomad fields AF and AC")


    args = parser.parse_args()
    return args


def process_vcf(args):
    # Load data
    # vcf_path = "/data/NICHD-core0/analysis/stratakis/poi-wes/workflows/variant-calling/results/annotated/ann.vcf.gz"
    # shortened_vcf_path = "results/annotated/subset.ann.vcf"
    # output_path = "/data/NICHD-core0/analysis/stratakis/poi-wes/workflows/variant-calling/results/analyzed"

    testing = False # Boolean to turn on if doing testing runs (runs with a subsetted vcf file)
    vcf = VCF(args.vcf, strict_gt=True)
    output_path = args.output
    if not os.path.exists(output_path):
        os.makedirs(output_path)


    print("Reading in data...")
    # W = Writer(output_path + '.vcf', vcf)
    samples = vcf.samples

    records = []
    for variant in vcf:
        ann_raw = variant.INFO.get("ANN")
        # These fields are assigned via SnpEff's documentation here: (https://pcingola.github.io/SnpEff/snpeff/inputoutput/)
        if ann_raw:
            ann_fields = ann_raw.split(",")[0].split("|")
            if ann_fields[2] not in ["MODERATE", "HIGH"]:
                continue
            allele = ann_fields[0] if len(ann_fields) > 0 else None
            var_type = ann_fields[1] if len(ann_fields) > 1 else None
            impact = ann_fields[2] if len(ann_fields) > 2 else None
            gene = ann_fields[3] if len(ann_fields) > 3 else None
            gene_id = ann_fields[4] if len(ann_fields) > 4 else None
            feature_type = ann_fields[5] if len(ann_fields) > 5 else None
            feature_id = ann_fields[6] if len(ann_fields) > 6 else None
            tscript_biotype = ann_fields[7] if len(ann_fields) > 7 else None
            ex_int_rank_total = ann_fields[8] if len(ann_fields) > 8 else None 
            HGVS_c = ann_fields[9] if len(ann_fields) > 9 else None
            HGVS_p = ann_fields[10] if len(ann_fields) > 10 else None
            cDNA_pos_rank = ann_fields[11] if len(ann_fields) > 11 else None
            CDS_pos_rank = ann_fields[12] if len(ann_fields) > 12 else None
            prot_pos_len = ann_fields[13] if len(ann_fields) > 13 else None
            feat_dist = ann_fields[14] if len(ann_fields) > 14 else None
            other = ann_fields[15] if len(ann_fields) > 15 else None

            record = {
                "CHROM": variant.CHROM,
                "POS": variant.POS,
                "REF": variant.REF,
                "ALT": variant.ALT,
                "DEPTH": variant.INFO.get("DP"),
                "GENOTYPES": variant.genotypes,
                "INFO_ANN": variant.INFO,
                "ANN": ann_raw,
                "ALLELE": allele,
                "VAR_TYPE": var_type,
                "IMPACT": impact,
                "GENE": gene,
                "GENE_ID": gene_id,
                "FEATURE_TYPE": feature_type,
                "FEATURE_ID": feature_id,
                "TSCRIPT_BIOTYPE": tscript_biotype,
                "EX_INT_RANK_TOTAL": ex_int_rank_total,
                "HGVS.C": HGVS_c,
                "HGVS.P": HGVS_p,
                "CDNA_POS_RANK": cDNA_pos_rank,
                "CDS_POS_RANK": CDS_pos_rank,
                "PROT_POS_LEN": prot_pos_len,
                "FEAT_DIST": feat_dist,
                "OTHER": other
            }

            # Add genotype info for each sample
            for i, sample in enumerate(samples):
                record[sample] = variant.gt_types[i]

            if args.gnomad:
                record["gnomad_AF"] = variant.INFO.get("AF_joint")
                record["gnomad_AC"] = variant.INFO.get("AC_joint")
                record["gnomad_AN"] = variant.INFO.get("AN_joint")
            records.append(record)


    vcf_df = pd.DataFrame(records)


    # Compute metrics for each gene and variant

    print("Computing variant metrics...")
    # Get list of unique genes in vcf
    genes = vcf_df['GENE'].unique()
    genes_df = pd.DataFrame(genes)
    genes_df.columns = ['gene']
    genes_df['num_vars'] = pd.Series([pd.NA] * len(genes_df), dtype='Int64')
    genes_df['sample_AF'] = np.NAN
    genes_df['sample_AC'] = pd.Series([pd.NA] * len(genes_df), dtype='Int64')
    genes_df['sample_AN'] = pd.Series([pd.NA] * len(genes_df), dtype='Int64')
    genes_df['vars'] = pd.Series([None] * len(genes_df), dtype='object')

    #TODO: Change this to an enumerate?
    index = 0
    for gene in genes:
        # Grab all rows/variants associated with the current gene
        cur_gene_df = vcf_df[vcf_df['GENE'] == gene]

        # Create list of dictionaries (each index corresponds to a variant for that gene, 
        # each index's dictionary will hold the info for that variant)
        vars = []
        total_alleles = 0
        total_alt_alleles = 0

        for _,var in cur_gene_df.iterrows():
            # In cyvcf2, the genotype numbers correspond as follows (even though their API says something else):
            # 0: homozygous reference
            # 1: heterozygous for alternate
            # 2: unknown (great ordering here by whoever made this up)
            # 3: homozygous for alternate
            var_alleles = 0
            # samples_gt = {sample: var[sample] for sample in samples}
            genotypes = var['GENOTYPES']
            samples_gt = {
                sample: "./." if genotypes[i][0] is None or genotypes[i][1] is None or genotypes[i][0] == -1 or genotypes[i][1] == -1
                    else f"{genotypes[i][0]}/{genotypes[i][1]}"
                for i, sample in enumerate(samples)
            }
            all_gt = var[samples] # Grab genotype info from all samples
            called_gt = all_gt[all_gt != 2] # Filter out 'unknown' variant calls
            var_samples = len(called_gt) # Number of samples with a genotype called for this variant
            var_alt_samples = len(called_gt[(called_gt==1) | (called_gt==3)]) # Number of samples with at least one called alternate allele for this variant
            var_alleles = len(called_gt) * 2 # Total alleles called for this variant
            alt_alleles = 0 # Total alternate alleles called for this variant
            for gt in called_gt:
                if gt == 1:
                    alt_alleles += 1
                elif gt == 3:
                    alt_alleles += 2
            total_alleles += var_alleles # Add to the total alleles called for the current gene
            total_alt_alleles += alt_alleles # Add to the total alt alleles called for the current gene
            var_AF = alt_alleles / var_alleles if var_alleles > 0 else 0 # Calculate the alternate allele frequency for this variant across all samples

            # Store info in list of variant info for current gene
            var_info = {"var_loc": f"{var['CHROM']}:{var['POS']}", "variant-id": f"{var['REF']}-{'/'.join(var['ALT'])}", "sample_AF": var_AF, "sample_AC": alt_alleles,"sample_AN": var_alleles,  "num_alt_samples": var_alt_samples, "total_called_samples": var_samples, "depth": var['DEPTH'], "var_impact": var['IMPACT'], "var_type": var["VAR_TYPE"], "genotypes": samples_gt}
            if args.gnomad:
                var_info["gnomad_AF"] = round(var['gnomad_AF'], 3)
                var_info["gnomad_AC"] = int(var['gnomad_AC']) if pd.notna(var['gnomad_AC']) else 0
                var_info["gnomad_AN"] = int(var['gnomad_AN']) if pd.notna(var['gnomad_AN']) else 0
            vars.append(var_info)

        # Compute alternate allele frequency across all samples
        gene_alt_freq = round(total_alt_alleles / total_alleles, 3) if total_alleles > 0 else 0

        genes_df.at[index, 'num_vars'] = len(vars)
        genes_df.at[index, 'sample_AF'] = gene_alt_freq
        genes_df.at[index, 'sample_AC'] = total_alt_alleles
        genes_df.at[index, 'sample_AN'] = total_alleles
        genes_df.at[index, 'vars'] = vars
        index += 1

    genes_df = genes_df.applymap(lambda x: float(f"{x:.4g}") if isinstance(x, float) and not pd.isnull(x) else x)

    # Expand variant column to create new files where each variant gets its own row
    vars_df = genes_df.drop(['num_vars', 'sample_AC','sample_AF', 'sample_AN'], axis=1)
    vars_df = vars_df.explode('vars')
    variant_details = pd.json_normalize(vars_df['vars'])
    vars_df= pd.concat([vars_df.drop(columns='vars').reset_index(drop=True), variant_details], axis=1)
    vars_df = vars_df.applymap(lambda x: float(f"{x:.4g}") if isinstance(x, float) and not pd.isnull(x) else x)
    vars_df.rename(columns={
        col: f"POF{col.split('POF')[-1]}_gt"
        for col in vars_df.columns
        if col.startswith('genotypes.POF')
    }, inplace=True)


    # Print the variants dataframe
    vars_df = vars_df.sort_values(by="sample_AC", ascending=False)
    with open(f'{output_path}/var_metrics.pkl', 'wb') as file:
        pickle.dump(vars_df, file)
    vars_df.to_csv(f'{output_path}/var_metrics.tsv', sep="\t", index=False)

    # Print the genes dataframe
    genes_final_df = genes_df.drop('vars', axis=1)
    genes_final_df = genes_final_df.sort_values(by="sample_AC", ascending=False)
    print("Saving to .pkl and .tsv...")
    if testing:
        with open(f'{output_path}/shortened_gene_var_metrics.pkl', 'wb') as file:
            pickle.dump(genes_final_df, file)

        genes_df.to_csv(f'{output_path}/shortened_gene_var_metrics.tsv', sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    else:
        with open(f'{output_path}/gene_var_metrics.pkl', 'wb') as file:
            pickle.dump(genes_final_df, file)

        genes_df.to_csv(f'{output_path}/gene_var_metrics.tsv', sep="\t", index=False)

    print("All done!")


if __name__ == "__main__":
    args = main()
    process_vcf(args)

