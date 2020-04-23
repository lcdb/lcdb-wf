#!/bin/bash

"""
This script take in a vcf output from four different somatic variant callers:
strelka, somaticsniper, mutect2, and varscan. It normalizes and left aligns the
calls from each caller, and splits out multiallelic calls. It also prepends the
caller name to the INFO and FORMAT fields for each variant caller.
"""


# Normalize and prepend FORMAT and INFO fields for Varscan VCF
(
bcftools norm 
-f snakemake.input[4]
--multiallelics -both
snakemake.input[0]
| sed "s/INFO\tFORMAT/VARSCAN_INFO\tVARSCAN_FORMAT/g" 
> snakemake.output[0]
) 

# Normalize and prepend FORMAT and INFO fields for Strelka VCF
(
bcftools norm 
-f snakemake.input[4]
--multiallelics -both
snakemake.input[1]
| sed "s/INFO\tFORMAT/STRELKA_INFO\tSTRELKA_FORMAT/g"
> snakemake.output[1]
)

#Normalize and prepend FORMAT and INFO fields for  SomaticSniper VCF
(
bcftools norm 
-f snakemake.input[4]
--multiallelics -both
snakemake.input[2]
| sed "s/INFO\tFORMAT/SOMATICSNIPER_INFO\tSOMATICSNIPER_FORMAT/g"
> snakemake.output[2]
)

#Normalize and prepend FORMAT and INFO fields for Mutect2 VCF
(
bcftools norm 
-f snakemake.input[4]
--multiallelics -both
snakemake.input[3]
| sed "s/INFO\tFORMAT/MUTECT2_INFO\tMUTECT2_FORMAT/g"
> snakemake.output[3]
)

