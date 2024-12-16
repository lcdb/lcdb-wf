import os
import pprint
from textwrap import dedent
from . import common
from . import multiqc
import yaml
import pandas as pd
from pandas.testing import assert_frame_equal

def test_config_loading(tmpdir):
    """
    Test the loading of configuration files, including inclusion and overriding of parameters
    from multiple YAML files and directories.

    Parameters
    ----------
    tmpdir : LocalPath
        Temporary directory provided by pytest for testing file-based workflows.

    Raises
    ------
    AssertionError
        If the loaded configuration does not match the expected dictionary.
    """
    f0 = tmpdir.mkdir('subdir').join('file0.yaml')
    dir_to_include = tmpdir.join('subdir')
    f0.write(dedent('''
    references:
      species_to_keep:
        tag_from_directory:
          fasta:
            url: "https://from_directory"

	# Will get overwritten by a specific file
        tag_from_file:
          fasta:
            url: "https://from_directory"

        # Will get overwritten by specific file, and then that will get
        # overwritten by the config
        tag_from_config:
          fasta:
            url: "https://from_directory"
    '''))
    f1 = tmpdir.join('subdir', 'file1.yaml')
    f1.write(dedent('''
    references:
      species2:
        tag_only_in_directory:
          fasta:
            url: ""
            indexes:
              - bowtie2
    '''))
    f2 = tmpdir.join('file1.yaml')
    f2.write(dedent('''
    references:
      species_to_keep:
        tag_from_file:
          fasta:
            url: "https://from_file"
        tag_from_config:
          fasta:
            url: "https://from_file"
    '''))
    f3 = tmpdir.join('file3.yaml')
    f3.write(dedent('''
    references_dir: "/data"
    references:
      species_to_keep:
        tag_from_config:
          fasta:
            url: "https://from_config"
    include_references:
      - {dir_to_include}
      - {f2}
    '''.format(dir_to_include=dir_to_include, f2=f2)))

    config = common.load_config(str(f3))

    assert config == {
        'references_dir': '/data',
        'include_references': [
            '{0}/subdir'.format(str(tmpdir)),
            '{0}/file1.yaml'.format(str(tmpdir)),
        ],
        'references': {
            'species_to_keep': {
                'tag_from_config': {
                    'fasta': {'url': 'https://from_config'}},
                'tag_from_directory': {
                    'fasta': {'url': 'https://from_directory'}},
                'tag_from_file': {
                    'fasta': {'url': 'https://from_file'}}
            },
            'species2': {
                'tag_only_in_directory': {
                    'fasta': {'indexes': ['bowtie2'], 'url': ''}}},
        },
    }


def test_multiqc_expr_module(tmpdir):
    """
    Test the MultiQC expression module for generating `expression_barchart_mqc.yaml`.

    Parameters
    ----------
    tmpdir : LocalPath
        Temporary directory provided by pytest for testing file-based workflows.

    Raises
    ------
    AssertionError
        If the generated YAML file does not match the expected structure and values.
    """
    known_truth_raw_counts = {'sample1': {'cold': 329, 'hop': 2977}, 'sample2': {'cold': 547,
                              'hop': 1847}, 'sample3': {'cold': 405, 'hop': 4024}, 'sample4':
                             {'cold': 290, 'hop': 3693}}
    known_truth_sizefactors = {'sample1': 0.898, 'sample2': 0.803,
                                'sample3': 1.164, 'sample4': 1.128}
    featurecounts_path = '../test/test_configs/test_featurecounts.txt'
    sampletable_path = '../test/test_configs/test_sampletable.tsv'
    processed_fc_path = tmpdir.mkdir('mqc_test').join('processed_featurecounts.tsv')
    sizefactors_path = tmpdir.join('sizefactors.tsv')
    yaml_path = tmpdir.join('expression_barchart_mqc.yaml')

    multiqc.process_featurecounts(featurecounts_path, sampletable_path, processed_fc_path)
    multiqc.estimate_sizefactors(processed_fc_path, sizefactors_path)
    counts_file_paths = populate_counts_files(known_truth_raw_counts, tmpdir)
    multiqc.make_regions_counts_yaml(sizefactors_path, counts_file_paths, known_truth_raw_counts['sample1'].keys(), yaml_path)
    test_yaml = yaml.safe_load(open(yaml_path, "r"))

    snakefile_yaml = {
      "data": {
          "sample1": {"cold": 366.437, "hop": 3315.754},
          "sample2": {"cold": 681.177, "hop": 2300.062},
          "sample3": {"cold": 348.072, "hop": 3458.376},
          "sample4": {"cold": 257.043, "hop": 3273.304},
      },
      "description": "The user-defined regions are: cold, hop",
      "id": "expression_barchart",
      "pconfig": {
          "id": "user_regions_counts_plot",
          "title": "Normalized counts in cold, hop",
          "xlab": "Samples",
          "ylab": "Normalized Counts",
      },
      "plot_type": "bargraph",
      "section_name": "Normalized Counts in User Defined Regions",
    }

    assert test_yaml == snakefile_yaml


def populate_counts_files(known_truth_raw_counts, tmpdir):
    """
    Populate temporary files with raw count data for each sample and region.

    Parameters
    ----------
    known_truth_raw_counts : dict
        Dictionary containing raw counts for each sample and region.
    tmpdir : LocalPath
        Temporary directory provided by pytest.

    Returns
    -------
    list of str
        List of file paths to the populated count files.
    """
    counts_file_paths = []
    for sample, regions in known_truth_raw_counts.items():
        for region, count in regions.items():
            subdir = tmpdir.join(sample)
            subdir.ensure(dir=True)
            file_path = subdir.join(f"{region}.counts.txt")
            file_path.write(str(count))
            counts_file_paths.append(str(file_path))
    return counts_file_paths


def test_process_featurecounts(tmpdir):
    """
    Test processing and normalization of feature counts using the MultiQC module.

    Parameters
    ----------
    tmpdir : LocalPath
        Temporary directory provided by pytest.

    Raises
    ------
    AssertionError
        If the processed and normalized feature counts do not match the expected values.
    """
    featurecounts_path = '../test/test_configs/test_featurecounts.txt'
    sampletable_path = '../test/test_configs/test_sampletable.tsv'
    processed_fc_path = tmpdir.mkdir('mqc_test').join('processed_featurecounts.tsv')
    sizefactors_path = tmpdir.join('sizefactors.tsv')
    normalized_path = tmpdir.join('normalized_featurecounts.tsv')
    multiqc.process_featurecounts(featurecounts_path, sampletable_path, processed_fc_path)
    multiqc.estimate_sizefactors(processed_fc_path, sizefactors_path)
    multiqc.normalize_featurecounts(processed_fc_path, sizefactors_path, normalized_path)
    test_normalized_fc = pd.read_csv(normalized_path, sep='\t')
    snakefile_normalized_fc = pd.read_csv('../test/test_configs/normalized_counts.tsv', sep='\t')
    assert_frame_equal(test_normalized_fc, snakefile_normalized_fc, atol=1e-3)
