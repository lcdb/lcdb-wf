import pandas as pd
import yaml
import os

# Read sample table
samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
units = pd.read_table(config["units"], dtype=str).set_index(["sample","unit"], drop=False)
units.index = units.index.set_levels([i for i in units.index.levels])

def preflight():
    """
    This helper function gets called at the top of the main Snakefile. It
    handles reading the config to see if references are provided externally, or
    if we are relying on lcdb-wf references. Returns variables containing
    filepaths of references to be used in rules. It will also perform some
    checks to make sure the config is not contradicting itself under certain
    configurations.
    """
    aln_index = []
    dbnsfp = []
    dictionary = []
    indexed = []
    known_sites = []
    reference = []

    # Handle reference names if LCDB-WF References is ran
    if config['ref']['use_references_workflow']:
        include: '../references/Snakefile'
        refdict = common.references_dict(config)
        reference = refdict[config['ref']['organism'][config['ref']['genome']['tag']]['genome']]
        aln = refdict[config['ref']['organism'][config['ref']['aligner']['tag']]['bwa']]
        aln_index = multiext(os.path.splitext(aln)[0], ".amb", ".ann", ".bwt", ".pac", ".sa")
        indexed = refdict[config['ref']['organism'][config['ref']['faidx']['tag']]['faidx']]
        if config['ref']['variation']['dbnsfp']:
            dbnsfp = refdict[config['ref']['organism']]['variation'][str(config['ref']['variation']['dbnsfp'] + '_' + config['ref']['genome']['build'])]
        else:
            dbnsfp = []
        if config['ref']['variation']['known']:
            known_sites = refdict[config['ref']['organism']][config['ref']['genome']['tag']][config['ref']['variation']['known']]
        else:
            known_sites = []
    else:
        known_sites = (
            config['ref']['paths']['known']
            if config['ref']['paths']['known']
            else []
        )
        reference = config['ref']['paths']['ref']
        indexed = (
            config['ref']['paths']['index']
            if config['ref']['paths']['index']
            else reference + '.fai'
        )
        dbnsfp = (
            config['ref']['paths']['dbnsfp']
            if config['ref']['paths']['dbnsfp']
            else []
        )
        aln_index = []

    # Handle dictionary name, stop the workflow if the fasta file is not named properly. Stop the workflow if there is no reference
    if reference == []:
        raise ValueError("You must supply a reference file to workflow.")
    if reference.endswith('.gz'):
        dictionary = '.'.join(reference.split('.')[:-2]) + '.dict'
    else:
        try:
            dictionary ='.'.join(reference.split('.')[:-1]) + '.dict'
            # If there is no exception, python will raise a TypeError trying to concatenate an empty list with str
        except TypeError:
            raise ValueError("There is something wrong with your reference extension. " 
                             "Please make sure your reference has an extension")
    # Stop the workflow easily if there is no known variation, but bqsr is set in the config
    if config['filtering']['bqsr'] == True:
        assert known_sites != [], 'Check your config.yaml. You are requiring that bqsr be run, but there is no known sites vcf'

    return aln_index, dbnsfp, dictionary, indexed, known_sites, reference


def get_contigs():
    """
    Helper function to read the contigs from the fasta index checkpoint rule.
    These contigs define the regions to split variant calling by for joint-calling.
    """
    with checkpoints.genome_index.get().output[0].open() as fai:
        ser = pd.read_table(fai, header=None, usecols=[0], dtype=str)
        ser = ser.squeeze()
        # TODO: make this less brittle, and better support non-Ensembl organisms
        # Remove all contigs that don't correspond to a chromosome
        ser = ser[ser.apply(lambda x: len(x)) <= 2]
        # Remove mitochondiral if specified in the config
        if config["processing"]["remove-mitochondrial"]:
            return ser[ser != "MT"]
        else:
            return ser


def get_fastq(wildcards):
    """
    Get fastq files of given sample-unit. Sample-unit structure is how technical replicates
    are handled. This is defined in the sampletable.
    """
    fastqs = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return "-R '@RG\\tID:{sample}\\tSM:{sample}\\tPL:{platform}'".format(
        sample=wildcards.sample,
        platform=units.loc[(wildcards.sample, wildcards.unit), "platform"],
    )


def get_recal_input(bai=False):
    """
    Handle providing bams the input of the bqsr rules.
    Read config options to determine the appropriate bam and bam index files.
    If we don't remove duplicates, return the sorted bams from map reads rule.
    If duplicates are removed, return the deduplicated bams from the mark duplicates rule.
    If a bed file is used in variant calling
    """
    # Case 1: no duplicate removal
    f = "results/mapped/{sample}-{unit}.sorted.bam"
    if config["processing"]["remove-duplicates"]:
        # Case 2: remove duplicates
        f = "results/dedup/{sample}-{unit}.bam"
    if bai:
        if config["processing"]["restrict-regions"]:
            # Case 3: need an index because random access is required
            f += ".bai"
            return f
        else:
            # Case 4: no index needed
            return []
    else:
        return f


def get_sample_bams(wildcards):
    """
    Get all aligned reads of given sample. Return the recal bams if bqsr is run
    otherwise return the dedup bams. We return all units for a given sample because
    we want to provide technical replicates to the variant calling rule where this is called
    """
    unitlist = units.loc[wildcards.sample].unit
    reslist = []
    if config['filtering']['bqsr']:
        reslist.extend(
            [
                "results/recal/{}-{}.bam".format(wildcards.sample, unit) for unit in unitlist
            ]
        )
    else:
        reslist.extend(
            [
                "results/dedup/{}-{}.bam".format(wildcards.sample, unit) for unit in unitlist
            ]
        )

    return reslist


def get_sample_unit_bams(wildcards):
    """
    Get all aligned reads of given sample. Unlike the function above, we return a single sample-unit combination per function call.
    This is because this function is used to QC rules like samtools-stats where we do not want to combine technical replicates.
    Return the recal bams if bqsr is run otherwise return the dedup bams
    """
    reslist = ''
    if config['filtering']['bqsr']:
        reslist = "results/recal/{sample}-{unit}.bam".format(sample=wildcards.sample, unit=wildcards.unit)
    else:
        reslist = "results/dedup/{sample}-{unit}.bam".format(sample=wildcards.sample, unit=wildcards.unit)
    return reslist


def get_regions_param(regions=config["processing"]["restrict-regions"], default=""):
    """
    If a captured regions bedfile is present, split the variant calling up into regions
    follwing GATK best practices
    """
    if regions:
        params = "--intervals '{}' ".format(regions)
        padding = config["processing"].get("region-padding")
        if padding:
            params += "--interval-padding {}".format(padding)
        return params
    return default


def get_call_variants_params(wildcards, input):
    """
    Calls the previous function to assemble the regions into interval lists
    along with any specified parameters for variant calling in the config
    """
    return (
        get_regions_param(
            regions=input.regions, default="--intervals {}".format(wildcards.contig)
        )
    )


def set_java_opts(resources):
    """
    Using the resources directive from the snakemake rule
    set the heap size. Request 75 percent of the requested
    mem_mb. The remaining 25 percent should be enough for 
    OS and other system processes that occur outside the shell command
    """
    heap = int(resources.mem_mb * 0.75)
    heap = int(heap / 1024)
    java_temp ='''"-Xmx{}g -Djava.io.tmpdir=$TMPDIR\"'''
    java_opts = java_temp.format(heap)
    return java_opts

def all_input_mutect2():
    """
    Format the input for the all rule for mutect2
    """
    comparisons = config['mutect2'].keys()
    return expand("results/mutect2_annotated_normed/{comp}.vcf.gz", comp=comparisons)


def names_for_somatic(wildcards):
    """
    Format the names into arguments to pass to mutect2.
    Mutect2 requires you to specify the names of the "normal" samples.
    There can be multiple normal samples in a single mutect2 call.
    Tumor samples do not need to be named. This will be done by reading
    from the config.
    """
    comp = wildcards.comp
    normals = config['mutect2'][comp]['normal']
    if not isinstance(normals, list):
        normals = [normals]
    return normals


def input_for_somatic(wildcards):
    """
    Format the bam input for mutect2 by reading from the config. 
    Technical replicates are separated and grouped. Returns a dictionary
    contains the reference genome, sequence dictionary, and input bams
    """
    comp = wildcards.comp
    normals = config['mutect2'][comp]['normal']
    if not isinstance(normals, list):
        normals = [normals]
    tumors = config['mutect2'][comp]['tumor']
    if not isinstance(tumors, list):
        tumors = [tumors]
    # Fill these lists with paths to tumor and normal files
    t_files = []
    n_files = []
    for i in range(len(tumors)):
        # Get the unit for each tumor sample
        unitlist = units.loc[tumors[i]].unit
        if config['filtering']['bqsr']:
            t_files.extend(
                [
                "results/recal/{}-{}.bam".format(tumors[i], unit) for unit in unitlist
                ]
            )
        else:
            t_files.extend(
                [
                "results/dedup/{}-{}.bam".format(tumors[i], unit) for unit in unitlist
                ]
            )
    # Do the same for Normals
    for i in range(len(normals)):
        unitlist = units.loc[normals[i]].unit
        if config['filtering']['bqsr']:
            n_files.extend(
                [
                "results/recal/{}-{}.bam".format(normals[i], unit) for unit in unitlist
                ]
            )
        else:
            n_files.extend(
                [
                "results/dedup/{}-{}.bam".format(normals[i], unit) for unit in unitlist
                ]
            )


    # Put all the input files needed into a dictionary to pass to the rule
    d = dict(
        ref=reference,
        normals=n_files,
        tumors=t_files,
        dict=dictionary,
        regions=(
            "results/called/{contig}.regions.bed".format(contig = wildcards.contig)
            if config["processing"]["restrict-regions"]
            else []
        ),
    )
    return d


def get_fai_nomenclature():
    """
    Helper function to get the nomenclature of the fasta index
    Returns True if the chr prefix is present, and False if it is absent
    """
    nom = False
    with checkpoints.genome_index.get().output[0].open() as fai:
        for line in fai:
            if line.startswith('chr'):
                nom = True
                break
    return nom


def get_bed_nomenclature(input):
    """
    Helper function to get the nomenclature of the bedfile
    Returns True if the chr prefix is present, and False if it is absent
    """
    nom = False
    with open(input.bed, 'r') as f:
        for line in f:
            if line.startswith('browser') or line.startswith('track'):
                continue
            if f.startswith('chr'):
                nom = True
                break
    return nom


# vim: ft=python
