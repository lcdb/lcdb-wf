import sys
import subprocess
top_level_dir = subprocess.run(["dirname $(dirname $(pwd))"], shell=True, capture_output=True, text=True).stdout.strip()
print("top level dir: ", top_level_dir)
sys.path.insert(0, top_level_dir)
import pytest
from textwrap import dedent
from lib import common, helpers, patterns_targets

# Make config object that can be re-used for any test
@pytest.fixture
def config(request):
    config_path = request.param
    config = common.load_config(config_path, test=True)
    return patterns_targets.RNASeqConfig(config, config.get('patterns', '../workflows/rnaseq/config/rnaseq_patterns.yaml'))

# Call helpers.detect_layout(), which implicitly tests common.is_paired_end()
# TODO: Make assertion condition NOT hard coded in to work with current example table
@pytest.mark.parametrize("config", ['../../workflows/rnaseq/config/config.yaml'], indirect=True)
def test_is_paired_end(config):
    is_paired = helpers.detect_layout(config.sampletable) == 'PE'
    assert not is_paired, f"Test failed, is_paired = {is_paired}"

def test_config_loading(tmpdir):
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

