import collections
import re
from itertools import product
import pandas as pd
from snakemake.shell import shell
from snakemake.io import expand, regex


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
    def update(d, u, c):
        for k, v in u.items():
            if isinstance(v, collections.Mapping):
                r = update(d.get(k, {}), v, c)
                d[k] = r
            else:
                if isinstance(fill, pd.DataFrame):
                    d[k] = list(set(expand(u[k], zip, **fill.to_dict('list'))))
                else:
                    d[k] = list(set(expand(u[k], c, **fill)))
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
    with open(scriptname, 'w') as fout:
        fout.write(string)
    if log:
        _log = '> {0} 2>&1'.format(log)
    else:
        _log = ""
    shell('Rscript {scriptname} {_log}')
