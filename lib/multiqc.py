def lcdbwf_samplename(x):
    """
    Processes sample names by removing specific substrings and extracting relevant parts.

    Args:
        x (list of str): A list of sample name strings.

    Returns:
        list of str: A list of processed sample names.
    """
    processed_names = []
    for name in x:
        # Remove unnecessary substrings and split to obtain the core sample name
        name = name.replace('data/rnaseq_samples/', '').replace('.cutadapt.bam', '')
        name_parts = name.split('/')
        processed_names.append(name_parts[0])
    return processed_names

def process_feature_counts(fc_path, sampletable_path, output_path, sample_func=lcdbwf_samplename, subset_counts=False):
    """
    Processes a feature counts table and matches it with a sample table, saving the processed table to the specified output path.

    Args:
        fc_path (str): Path to the feature counts file.
        sampletable_path (str): Path to the sample table TSV file.
        output_path (str): Path to save the processed counts DataFrame.
        sample_func (function): Function to process sample names from counts data.
        subset_counts (bool): If True, subsets counts to samples present in the sample table.

    Returns:
        None
    """

    import pandas as pd

    # Load counts, set 'Geneid' as index, remove metadata columns, and process sample names
    m = pd.read_csv(fc_path[0], sep='\t', comment='#').set_index('Geneid').iloc[:, 5:]
    x = sample_func(m.columns.tolist())
    sampletable = pd.read_csv(sampletable_path, sep='\t')
    samplenames = sampletable.iloc[:, 0].tolist()

    # Identify samples mismatched between counts and sampletable
    counts_not_sampletable = set(x) - set(samplenames)
    sampletable_not_counts = set(samplenames) - set(x)

    # Handle mismatches and subset counts if specified
    if counts_not_sampletable:
        if not subset_counts:
            raise ValueError(
                f"The following samples are in the counts data but not in the sample table. "
                f"Consider using subset_counts=True to remove them from the counts: "
                f"{', '.join(counts_not_sampletable)}"
            )
        else:
            indices_to_keep = [i for i, name in enumerate(x) if name in samplenames]
            x = [name for name in x if name in samplenames]
            m = m.iloc[:, indices_to_keep]

    if sampletable_not_counts:
        raise ValueError(
            f"The following samples are in the sample table but not in the counts data. "
            f"Check sample_func? {', '.join(sampletable_not_counts)}"
        )

    # Align columns in m to sampletable order and save to output path
    m.columns = x
    m = m[samplenames]
    m.to_csv(output_path[0], sep='\t')

def geometric_mean(arr):
    """
    Calculates the geometric mean of non-zero elements of an array.

    Args:
        arr (array-like): Array of values.

    Returns:
        float: Geometric mean of non-zero elements in the array.
    """
    import numpy as np
    return np.exp(np.mean(np.log(arr)))

def estimate_size_factors(counts, size_factors):
    """
    Estimates size factors for normalization using the median of ratios of counts to geometric means.

    Args:
        counts (str): Path to the counts TSV file.
        size_factors (str): Path to save the estimated size factors as a TSV file.

    Returns:
        None
    """

    import pandas as pd
    import numpy as np

    # Load counts, remove rows with zeros, and compute geometric means for each gene
    counts_df = pd.read_csv(counts[0], sep='\t', index_col=0)
    counts_df = counts_df[counts_df > 0].dropna()
    geometric_means = counts_df.apply(geometric_mean, axis=1)

    # Compute size factors by taking median of ratios of counts to geometric means
    ratios = counts_df.div(geometric_means, axis=0)
    size_factors_df = ratios.apply(np.median, axis=0)

    # Format size factors as DataFrame, reorder columns, and save
    size_factors_df = pd.DataFrame(size_factors_df, columns=['sizeFactors'])
    size_factors_df['samplename'] = size_factors_df.index
    size_factors_df = size_factors_df[['samplename', 'sizeFactors']].reset_index(drop=True)
    size_factors_df.to_csv(size_factors[0], sep='\t', header=True, index=False)

def regions_yaml(size_factors, counts, regions, yaml_file):
    """
    Generates a YAML file with normalized counts for specified regions for visualization in MultiQC.

    Args:
        size_factors (str): Path to the size factors TSV file.
        counts (list of str): List of paths to counts files for different regions.
        yaml (str): Path to save the generated YAML file.

    Returns:
        None
    """

    import yaml
    import pandas as pd
    import os
    from collections import defaultdict

    # Load size factors and initialize dictionary for normalized counts
    size_factors_df = pd.read_csv(size_factors[0], sep='\t', index_col='samplename')
    normalized_counts = defaultdict(dict)

    # Calculate normalized counts for each region and store in dictionary
    for f in counts:
        parts = f.strip().split(os.sep)
        sample = parts[2]
        region = parts[3].replace('.counts.txt', '')

        with open(f) as infile:
            count = int(infile.read())

        size_factor = size_factors_df.loc[sample, 'sizeFactors']
        normalized_counts[sample][region] = float(count) / float(size_factor)

    # Prepare YAML structure with plot metadata for MultiQC and write to file
    data = {sample: regions for sample, regions in normalized_counts.items()}
    y = {
        'id': 'expression_barchart',
        'section_name': 'Normalized Counts in User Defined Regions',
        'description': f"The user-defined regions are: {', '.join(regions[0])}",
        'plot_type': 'bargraph',
        'pconfig': {
            'id': 'user_regions_counts_plot',
            'title': f"Normalized counts in {', '.join(regions[0])}",
            'ylab': 'Normalized Counts',
            'xlab': 'Samples'
        },
        'data': data
    }

    with open(yaml_file[0], 'w') as outfile:
        yaml.dump(y, outfile, default_flow_style=False)

def normalize_counts(counts_path, size_factors_path, output_path):
    """
    Normalizes the counts table by size factors for each sample.

    Args:
        counts_path (str): Path to the counts TSV file.
        size_factors_path (str): Path to the size factors TSV file.
        output_path (str): Path to save the normalized counts TSV file.

    Returns:
        None
    """

    import pandas as pd

    # Read counts and size factors
    counts_df = pd.read_csv(counts_path[0], sep='\t', index_col=0)
    size_factors_df = pd.read_csv(size_factors_path[0], sep='\t', index_col='samplename')

    # Normalize each sample's counts by its size factor
    for sample in counts_df.columns:
        if sample in size_factors_df.index:
            counts_df[sample] = counts_df[sample] / size_factors_df.loc[sample, 'sizeFactors']
        else:
            raise ValueError(f"Size factor for sample '{sample}' not found in size factors file.")

    # Save the normalized counts table
    counts_df.to_csv(output_path[0], sep='\t')

