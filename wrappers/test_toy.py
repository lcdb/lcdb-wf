import os
from textwrap import dedent
import pytest
import utils

# Each module has a config dict
config = dict()


def generic_fixture(key, mapping, factory):
    """
    Tries to handle as much of the magic as possible.

    Parameters
    ----------
    key : str
        Key into the module-level config dict

    mapping : dict
        Maps paths from fixtures to input files expected by the snakefile

    tmpdir : str
        Path to temporary dir, usually created by utils.tmpdir_for_func

    Returns
    -------
    After a successful Snakemake run, returns the dictionary of the config's
    `output` key but with paths fixed to be relative to tmpdir. This returned
    dict is ready to be used as a fixture by test functions.
    """
    conf = config[key]
    tmpdir = utils.tmpdir_for_func(factory)
    input_data_func = utils.symlink_in_tempdir(mapping)
    utils.run(utils.dpath(conf['wrapper']), conf['snakefile'], None, input_data_func, tmpdir)
    output = conf['output'].copy()
    for k, v in output.items():
        output[k] = os.path.join(tmpdir, v)
    return output


# In order for the doc generation to find this config info without re-running
# all tests, it needs to be in the module-level dict. It similarly can't be
# added during the fixture function's runtime.
#
# However, the mapping and tmpdir must be provided by the function, so the
# config and the function are tightly coupled.
#
# So we add the item to the dictionary here, right above the function that will
# be using it to keep them tightly coupled in the file.
config['hisat2_index'] = dict(
    description="Basic example of generating a hisat2 index",
    wrapper="../wrappers/hisat2/build",
    snakefile="""
      rule hisat2_build:
          input:
              fasta="2L.fa"
          output:
              index=expand("hisat2_index/assembly.{n}.ht2", n=range(1,9))
          log: "hisat.log"
          wrapper: "file://wrapper"
    """,
    output={'prefix': 'hisat2_index/assembly'}
)


# All the hard work is done in the config and in generic_fixture(). Now we just
# need to set up the correct mapping of fixtures to input files.
@pytest.fixture(scope='module')
def hisat2_index(tmpdir_factory, dm6_fa):
    mapping = {dm6_fa: '2L.fa'}
    return generic_fixture('hisat2_index', mapping, tmpdir_factory)

# The actual test.
def test_index(hisat2_index):
    assert os.path.exists(hisat2_index['prefix'] + '.1.ht2')


def extract_examples_for_wrapper(wrapper):
    """
    Returns the examples for the wrapper in markdown format.

    Parameters
    ----------
    wrapper : str
        Expected to be the value of one of the config dict's `wrapper` keys.
    """
    markdown = []
    for k, v in config.items():
        if v['wrapper'] != wrapper:
            continue
        snakefile = dedent(v['snakefile'])
        markdown.append(
            dedent(
                """
                {}

                ```python""".format(v['description'])))
        markdown.append(snakefile)
        markdown.append("```")
    return "\n".join(markdown)
