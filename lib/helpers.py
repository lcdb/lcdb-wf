import collections
import re
from itertools import product
import pandas as pd
from snakemake.shell import shell
from snakemake.io import expand, regex
from lib import common


class ConfigurationError(Exception):
    pass


def detect_layout(sampletable):
    """
    Identifies whether a sampletable represents single-end or paired-end reads.

    Raises NotImplementedError if there's a mixture.
    """
    is_pe = [common.is_paired_end(sampletable, s) for s in sampletable.iloc[:, 0]]
    if all(is_pe):
        return "PE"
    elif not any(is_pe):
        return "SE"
    else:
        p = sampletable.iloc[is_pe, 0].to_list()
        s = sampletable.iloc[[not i for i in is_pe], 0].to_list()
        if len(p) > len(s):
            report = f"SE samples: {s}"
        else:
            report = f"PE samples: {p}"
        raise ValueError(f"Only a single layout (SE or PE) is supported. {report}")


def fill_patterns(patterns, fill, combination=product):
    """
    Fills in a dictionary of patterns with the dictionary or DataFrame `fill`.

    >>> patterns = dict(a='{sample}_R{N}.fastq')
    >>> fill = dict(sample=['one', 'two'], N=[1, 2])
    >>> sorted(fill_patterns(patterns, fill)['a'])
    ['one_R1.fastq', 'one_R2.fastq', 'two_R1.fastq', 'two_R2.fastq']

    >>> patterns = dict(a='{sample}_R{N}.fastq')
    >>> fill = dict(sample=['one', 'two'], N=[1, 2])
    >>> sorted(fill_patterns(patterns, fill, zip)['a'])
    ['one_R1.fastq', 'two_R2.fastq']

    >>> patterns = dict(a='{sample}_R{N}.fastq')
    >>> fill = pd.DataFrame({'sample': ['one', 'two'], 'N': [1, 2]})
    >>> sorted(fill_patterns(patterns, fill)['a'])
    ['one_R1.fastq', 'two_R2.fastq']

    """
    # In recent Snakemake versions (e.g., this happens in 5.4.5) file patterns
    # with no wildcards in them are removed from expand when `zip` is used as
    # the combination function.
    #
    # For example, in 5.4.5:
    #
    #   expand('x', zip, d=[1,2,3]) == []
    #
    # But in 4.4.0:
    #
    #   expand('x', zip, d=[1,2,3]) == ['x', 'x', 'x']

    def update(d, u, c):
        for k, v in u.items():
            if isinstance(v, collections.abc.Mapping):
                r = update(d.get(k, {}), v, c)
                d[k] = r
            else:
                if isinstance(fill, pd.DataFrame):
                    d[k] = list(set(expand(u[k], zip, **fill.to_dict("list"))))
                else:
                    d[k] = list(set(expand(u[k], c, **fill)))
            if not d[k]:
                d[k] = [u[k]]
        return d

    d = {}
    return update(d, patterns, combination)


def extract_wildcards(pattern, target):
    """
    Return a dictionary of wildcards and values identified from `target`.

    Returns None if the regex match failed.

    Parameters
    ----------
    pattern : str
        Snakemake-style filename pattern, e.g. ``{output}/{sample}.bam``.

    target : str
        Filename from which to extract wildcards, e.g., ``data/a.bam``.

    Examples
    --------
    >>> pattern = '{output}/{sample}.bam'
    >>> target = 'data/a.bam'
    >>> expected = {'output': 'data', 'sample': 'a'}
    >>> assert extract_wildcards(pattern, target) == expected
    >>> assert extract_wildcards(pattern, 'asdf') is None
    """
    m = re.compile(regex(pattern)).match(target)
    if m:
        return m.groupdict()


def rscript(string, scriptname, log=None):
    """
    Saves the string as `scriptname` and then runs it

    Parameters
    ----------
    string : str
        Filled-in template to be written as R script

    scriptname : str
        File to save script to

    log : str
        File to redirect stdout and stderr to. If None, no redirection occurs.
    """
    with open(scriptname, "w") as fout:
        fout.write(string)
    if log:
        _log = "> {0} 2>&1".format(log)
    else:
        _log = ""
    shell("Rscript {scriptname} {_log}")


def check_unique_fn(df):
    """
    Raises an error if the fastq filenames are not unique
    """
    fns = df["orig_filename"]
    if "orig_filename_R2" in df.columns:
        fns = pd.concat([fns, df["orig_filename_R2"]])
    if len(fns.unique()) < len(fns):
        raise ValueError("Fastq filenames non unique, check the sampletable\n")


def check_unique_samplename(df):
    """
    Raises an error if the samplenames are not unique
    """
    ns = df.index
    if len(ns.unique()) < len(ns):
        raise ConfigurationError("Samplenames non unique, check the sampletable\n")


def preflight(config):
    """
    Performs verifications on config and sampletable files

    Parameters
    ----------
    config: yaml config object
    """
    sampletable = pd.read_table(config["sampletable"], index_col=0, comment="#")
    check_unique_samplename(sampletable)
    if "orig_filename" in sampletable.columns:
        check_unique_fn(sampletable)


def rnaseq_preflight(c):
    if "kallisto" not in c.config:
        raise ConfigurationError(
            """
            Starting in v1.8, an additional 'kallisto' argument is expected
            in the config file. Note that in the future this may be
            automatically included, but for now please add the following to the
            config, where 'tagname' is the tag for the reference of interest:

            kallisto:
              tag: "tagname"
            """
        )


def chipseq_preflight(c):
    pass


def strand_arg_lookup(config, lookup):
    """
    Given a config object and lookup dictionary, confirm that the config has
    correctly specified strandedness and then return the value for that key.
    """
    if not config.stranded:
        raise ConfigurationError(
            "Starting in v1.8, 'stranded' is required in the config file. "
            "Values can be 'unstranded', 'fr-firststrand' (R1 aligns antisense to original transcript), "
            "or 'fr-secondstrand' (R1 aligns sense to original transcript). If you are not sure, "
            "run the workflow with only the 'strand_check' rule, like "
            "'snakemake -j 5 strand_check'."
        )
    if config.stranded not in lookup:
        keys = list(lookup.keys())
        raise KeyError(f"'{config.stranded}' not one of {keys}")
    return lookup[config.stranded]
