import binascii
import collections
import contextlib
import gzip
import os
import re
import subprocess
import warnings
from collections.abc import Iterable
from itertools import product

import pandas
import pandas as pd
import yaml
from Bio import SeqIO
from snakemake.io import expand, regex_from_filepattern
from snakemake.shell import shell

# Small helper functions


def resolve_name(name):
    """
    Imports a specific object from a dotted path and returns just that object.

    From nose.utils.resolve_name (with the logging parts taken out) which in
    turn is from unittest.TestLoader.loadTestByName
    """
    parts = name.split(".")
    parts_copy = parts[:]
    while parts_copy:
        try:
            module = __import__(".".join(parts_copy))
            break
        except ImportError:
            del parts_copy[-1]
            if not parts_copy:
                raise
    parts = parts[1:]
    obj = module
    for part in parts:
        obj = getattr(obj, part)
    return obj


@contextlib.contextmanager
def temp_env(env):
    """
    Context manager to temporarily set os.environ.
    """
    env = dict(env)
    orig = os.environ.copy()
    _env = {k: str(v) for k, v in env.items()}
    os.environ.update(_env)
    try:
        yield
    finally:
        os.environ.clear()
        os.environ.update(orig)


def flatten(iter, unlist=False):
    """
    Flatten an arbitrarily nested iterable whose innermost items are strings
    into a flat list of strings.

    Parameters
    ----------
    iter : iterable

    unlist : bool
        If True, convert single-item lists into a bare string
    """
    if isinstance(iter, dict):
        iter = iter.values()

    def gen():
        for item in iter:
            if isinstance(item, dict):
                item = item.values()
            if isinstance(item, Iterable) and not isinstance(item, str):
                yield from flatten(item)
            else:
                yield item

    results = list(gen())
    if unlist and len(results) == 1:
        return results[0]
    return results


def test_flatten():
    assert sorted(
        flatten(
            {
                "a": {
                    "b": {
                        "c": ["a", "b", "c"],
                    },
                },
                "x": ["e", "f", "g"],
                "y": {"z": "d"},
            }
        )
    ) == ["a", "b", "c", "d", "e", "f", "g"]

    assert flatten("a", True) == "a"
    assert flatten(["a"], True) == "a"
    assert flatten("a") == ["a"]
    assert flatten(["a"]) == ["a"]


def updatecopy(orig, update_with, keys=None, override=False):
    """
    Update a copy of a dictionary, with a bit more control than the built-in
    dict.update.

    Parameters
    -----------

    orig : dict
        Dict to update

    update_with : dict
        Dict with new values

    keys : list or None
        If not None, then only consider these keys in `update_with`. Otherwise
        consider all.

    override : bool
        If True, then this is similar to `dict.update`, except only those keys
        in `keys` will be considered. If False (default), then if a key exists
        in both `orig` and `update_with`, no updating will occur so `orig` will
        retain its original value.
    """
    d = orig.copy()
    if keys is None:
        keys = update_with.keys()
    for k in keys:
        if k in update_with:
            if k in d and not override:
                continue
            d[k] = update_with[k]
    return d


def update_recursive(orig, update_with):
    """
    Recursively update one dict with another.

    From https://stackoverflow.com/a/3233356

    >>> orig = {'a': {'b': 1, 'c': 2, 'd': [7, 8, 9]}}
    >>> update_with = {'a': {'b': 5}}
    >>> expected = {'a': {'b': 5, 'c': 2, 'd': [7, 8, 9]}}
    >>> result = update_recursive(orig, update_with)
    >>> assert result == expected, result

    >>> update_with = {'a': {'d': 1}}
    >>> result = update_recursive(orig, update_with)
    >>> expected = {'a': {'b': 5, 'c': 2, 'd': 1}}
    >>> result = update_recursive(orig, update_with)
    >>> assert result == expected, result
    """
    for k, v in update_with.items():
        if isinstance(v, collections.abc.Mapping):
            orig[k] = update_recursive(orig.get(k, {}), v)
        else:
            orig[k] = v
    return orig


def boolean_labels(names, idx, mapping={True: "AND", False: "NOT"}, strip="AND_"):
    """
    Creates labels for boolean lists.

    For example:

    >>> names = ['exp1', 'exp2', 'exp3']
    >>> idx = [True, True, False]
    >>> boolean_labels(names, idx)
    'exp1_AND_exp2_NOT_exp3'

    Parameters
    ----------

    names : list
        List of names to include in output

    idx : list
        List of booleans, same size as `names`

    mapping : dict
        Linking words to use for True and False

    strip : str
        Strip this text off the beginning of labels.

    given a list of names and a same-size boolean, return strings like

    a_NOT_b_AND_c

    or

    a_AND_b_AND_c_NOT_d_AND_e
    """
    s = []
    for n, x in zip(names, idx):
        s.append(mapping[x] + "_" + n)
    s = "_".join(s)
    if s.startswith(strip):
        s = s.replace(strip, "", 1)
    return s


def make_relative_symlink(target, linkname):
    """
    Helper function to create a relative symlink.

    Changes to the dirname of the linkname and figures out the relative path to
    the target before creating the symlink.
    """
    linkdir = os.path.dirname(linkname)
    relative_target = os.path.relpath(target, start=linkdir)
    linkbase = os.path.basename(linkname)
    if not os.path.exists(linkdir):
        shell("mkdir -p {linkdir}")
    shell(f"cd {linkdir}; ln -sf {relative_target} {linkbase}")


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
    m = re.compile(regex_from_filepattern(pattern)).match(target)
    if m:
        return m.groupdict()


def _is_gzipped(fn):
    """
    Filename-independent method of checking if a file is gzipped or not. Uses
    the magic number.

    xref https://stackoverflow.com/a/47080739
    """
    with open(fn, "rb") as f:
        return binascii.hexlify(f.read(2)) == b"1f8b"


def openfile(tmp, mode):
    """
    Returns an open file handle; auto-detects gzipped files.
    """
    if _is_gzipped(tmp):
        return gzip.open(tmp, mode)
    else:
        return open(tmp, mode)


def gzipped(tmpfiles, outfile):
    """
    Cat-and-gzip a list of uncompressed files into a compressed output file.
    """
    with gzip.open(outfile, "wt") as fout:
        for f in tmpfiles:
            with open(f) as infile:
                for line in infile:
                    fout.write(line)


def cat(tmpfiles, outfile):
    """
    Simple concatenation of files.

    Note that gzipped files can be concatenated as-is without un- and re-
    compressing.
    """
    shell(f"cat {tmpfiles} > {outfile}")


def is_paired_end(sampletable, sample):
    """
    Inspects the sampletable to see if the sample is paired-end or not

    Parameters
    ----------
    sampletable : pandas.DataFrame
        Contains a "layout" or "LibraryLayout" column (but not both). If the
        lowercase value is "pe" or "paired", consider the sample paired-end.
        Otherwise consider single-end.

    sample : str
        Assumed to be found in the first column of `sampletable`
    """
    # We can't fall back to detecting PE based on two fastq files provided for
    # each sample when it's an SRA sampletable (which only has SRR accessions).
    #
    # So detect first detect if SRA sampletable based on presence of "Run"
    # column and all values of that column starting with "SRR", and then raise
    # an error if the Layout column does not exist.

    if "Run" in sampletable.columns:
        if all(sampletable["Run"].str.startswith("SRR")):
            if (
                "Layout" not in sampletable.columns
                and "layout" not in sampletable.columns
            ):
                raise ValueError(
                    "Sampletable appears to be SRA, but no 'Layout' column "
                    "found. This is required to specify single- or paired-end "
                    "libraries."
                )

    row = sampletable.set_index(sampletable.columns[0]).loc[sample]
    if "orig_filename_R2" in row:
        return True
    if "layout" in row and "LibraryLayout" in row:
        raise ValueError("Expecting column 'layout' or 'LibraryLayout', " "not both")
    try:
        return row["layout"].lower() in ["pe", "paired"]
    except KeyError:
        pass
    try:
        return row["LibraryLayout"].lower() in ["pe", "paired"]
    except KeyError:
        pass
    return False


def fill_r1_r2(sampletable, pattern, r1_only=False):
    """
    Returns a function intended to be used as a rule's input function.

    The returned function, when provided with wildcards, will return one or two
    rendered versions of a pattern depending on SE or PE respectively.
    Specifically, given a pattern (which is expected to contain a placeholder
    for "{sample}" and "{n}"), look up in the sampletable whether or not it is
    paired-end.

    Parameters
    ----------

    sampletable : pandas.DataFrame
        Contains a "layout" column with either "SE" or "PE", or "LibraryLayout"
        column with "SINGLE" or "PAIRED". If column does not exist, assume SE.

    pattern : str
        Must contain at least a "{sample}" placeholder.

    r1_only : bool
        If True, then only return the file for R1 even if PE is configured.
    """

    def func(wc):
        try:
            wc.sample
        except AttributeError:
            raise ValueError(
                'Need "{{sample}}" in pattern ' '"{pattern}"'.format(pattern=pattern)
            )
        n = [1]
        if is_paired_end(sampletable, wc.sample) and not r1_only:
            n = [1, 2]
        res = expand(pattern, sample=wc.sample, n=n)
        return res

    return func


def pluck(obj, kv):
    """
    For a given dict or list that somewhere contains keys `kv`, return the
    values of those keys.

    Named after the dplyr::pluck, and implemented based on
    https://stackoverflow.com/a/1987195
    """
    if isinstance(obj, list):
        for i in obj:
            for x in pluck(i, kv):
                yield x
    elif isinstance(obj, dict):
        if kv in obj:
            yield obj[kv]
        for j in obj.values():
            for x in pluck(j, kv):
                yield x


# Functions for conveniently working with resources


def autobump(*args, **kwargs):
    """
    Used to automatically bump resources depending on how many times the job
    was attempted. This will return a function that is appropriate to use for
    an entry in Snakemake's `resources:` directive::

        rule example:
            input: "a.txt"
            resources:
                mem_mb=autobump(gb=10),
                runtime=autobump(hours=2, increment_hours=10)

    Values can be specified in multiple ways.

    A single number will be provided as the resource, and will be used to
    increment each time. For example, this is the equivalent of 10 GB for the
    first attempt, and 20 GB for the second:

    >>> f = autobump(1024 * 10)
    >>> f(None, 1)
    10240

    Adding a second unnamed argument will use it as a value to increment by for
    each subsequent attempt. This will use 10 GB for the first attempt, and 110
    GB for the second attempt.

    >>> f = autobump(1024 * 10, 1024 * 100)
    >>> f(None, 1)
    10240

    >>> f(None, 2)
    112640

    Instead of bare numbers, keyword arguments can be used for more convenient
    specification of units. The above two examples can also take this form:

    >>> f = autobump(gb=10)
    >>> f(None, 1)
    10240

    >>> f = autobump(gb=10, increment_gb=100)
    >>> f(None, 2)
    112640


    Units can be minutes, hours, days, mb, gb, or tb. For example:

    >>> f = autobump(hours=2, increment_hours=5)
    >>> f(None, 2)
    420

    """
    multiplier = {
        "mb": 1,
        "minutes": 1,
        "gb": 1024,
        "hours": 60,
        "days": 1440,
        "tb": 1024 * 1024,
    }
    units = list(multiplier.keys())

    if args and kwargs:
        raise ValueError(
            "Mixture of unnamed and keyword arguments not supported with autobump()"
        )

    if len(kwargs) > 2:
        raise ValueError("Only 2 kwargs allowed for autobump()")

    elif len(args) == 1 and not kwargs:
        baseline_converted = args[0]
        increment_converted = baseline_converted

    elif len(args) == 2 and not kwargs:
        baseline_converted, increment_converted = args

    elif len(kwargs) <= 2:
        baseline_kwargs = [k for k in kwargs.keys() if k in units]
        if len(baseline_kwargs) != 1:
            raise ValueError(
                "Multiple baseline kwargs found. Do you need to change one to have an 'increment_' prefix?"
            )

        baseline_kwarg = baseline_kwargs[0]
        baseline_value = kwargs[baseline_kwarg]
        baseline_unit = baseline_kwarg

        increment_kwargs = [k for k in kwargs if k.startswith("increment_")]
        if increment_kwargs:
            assert len(increment_kwargs) == 1
            increment_kwarg = increment_kwargs[0]
            increment_value = kwargs[increment_kwarg]
            increment_unit = increment_kwarg.split("_")[-1]
        else:
            increment_value = baseline_value
            increment_unit = baseline_unit

        if baseline_unit not in multiplier:
            raise ValueError(
                f"Baseline unit {baseline_unit} not in valid units {units}"
            )
        if increment_unit not in multiplier:
            raise ValueError(
                f"Increment unit {increment_unit} not in valid units {units}"
            )

        baseline_converted = baseline_value * multiplier[baseline_unit]
        increment_converted = increment_value * multiplier[increment_unit]

    else:
        raise ValueError(f"Unhandled args and kwargs: {args}, {kwargs}")

    def f(wildcards, attempt):
        return baseline_converted + (attempt - 1) * increment_converted

    return f


def gb(size_in_gb):
    return 1024 * size_in_gb


def hours(time_in_hours):
    return time_in_hours * 60


# Config parsing and handling


class ConfigurationError(Exception):
    pass


def detect_layout(sampletable):
    """
    Identifies whether a sampletable represents single-end or paired-end reads.

    Raises NotImplementedError if there's a mixture.
    """
    is_pe = [is_paired_end(sampletable, s) for s in sampletable.iloc[:, 0]]
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
    Fills in a dictionary of patterns with the dictionary `fill`.

    >>> patterns = dict(a='{sample}_R{N}.fastq')
    >>> fill = dict(sample=['one', 'two', 'three'], N=[1, 2])
    >>> sorted(fill_patterns(patterns, fill)['a'])
    ['one_R1.fastq', 'one_R2.fastq', 'three_R1.fastq', 'three_R2.fastq', 'two_R1.fastq', 'two_R2.fastq']

    If using `zip` as a combination, checks to ensure all values in `fill` are
    the same length to avoid truncated output.

    This fails:

    >>> patterns = dict(a='{sample}_R{N}.fastq')
    >>> fill = dict(sample=['one', 'two', 'three'], N=[1, 2])
    >>> sorted(fill_patterns(patterns, fill, zip)['a']) # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    ...
    ValueError: {'sample': ['one', 'two', 'three'], 'N': [1, 2]} does not have the same number of entries for each key

    But this works:

    >>> patterns = dict(a='{sample}_R{N}.fastq')
    >>> fill = dict(sample=['one', 'one', 'two', 'two', 'three', 'three'], N=[1, 2, 1, 2, 1, 2])
    >>> sorted(fill_patterns(patterns, fill, zip)['a'])
    ['one_R1.fastq', 'one_R2.fastq', 'three_R1.fastq', 'three_R2.fastq', 'two_R1.fastq', 'two_R2.fastq']

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

    if combination == zip:
        lengths = set([len(v) for v in fill.values()])
        if len(lengths) != 1:
            raise ValueError(
                f"{fill} does not have the same number of entries for each key"
            )

    def update(d, u, c):
        for k, v in u.items():
            if isinstance(v, collections.abc.Mapping):
                r = update(d.get(k, {}), v, c)
                d[k] = r
            else:  # not a dictionary, so we're at a leaf
                if isinstance(fill, pd.DataFrame):
                    d[k] = list(set(expand(u[k], zip, **fill.to_dict("list"))))
                else:
                    d[k] = list(set(expand(u[k], c, **fill)))
            if not d[k]:
                d[k] = [u[k]]
        return d

    d = {}
    return update(d, patterns, combination)


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
    pass


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


def filter_fastas(tmpfiles, outfile, pattern):
    """
    Extract records from fasta file(s) given a search pattern.

    Given input gzipped FASTAs, create a new gzipped fasta containing only
    records whose description matches `pattern`.

    Parameters
    ----------
    tmpfiles : list
        gzipped fasta files to look through

    outfile : str
        gzipped output fastq file

    pattern : str
        Look for this string in each record's description

    """

    def gen():
        for tmp in tmpfiles:
            handle = gzip.open(tmp, "rt")
            parser = SeqIO.parse(handle, "fasta")
            for rec in parser:
                if pattern not in rec.description:
                    continue
                rec.seq = rec.seq.back_transcribe()
                rec.description = rec.name
                yield rec

    with gzip.open(outfile, "wt") as fout:
        SeqIO.write(gen(), fout, "fasta")


def twobit_to_fasta(tmpfiles, outfile):
    """
    Converts .2bit files to fasta.

    Parameters
    ----------
    tmpfiles : list
        2bit files to convert

    outfile : str
        gzipped output fastq file
    """
    # Note that twoBitToFa doesn't support multiple input files, but we want to
    # support them with this function
    lookup = {i: i + ".fa" for i in tmpfiles}
    for i in tmpfiles:
        fn = lookup[i]
        shell("twoBitToFa {i} {fn}")

    # Make sure we retain the order of the originally-provided files from the
    # config when concatenating.
    fastas = [lookup[i] for i in tmpfiles]
    shell("cat {fastas} | gzip -c > {outfile}")
    shell("rm {fastas}")


def download_and_postprocess(urls, postprocess, outfile, log):
    """
    Many reference files cannot be used as-is and need to be modified.

    This function supports providing one or more URLs, and any postprocess
    functions to get the reference files usable.

    Parameters
    ----------
    urls : str or list
        URL(s) to download. Can be a list, in which case they will be concatenated.

    postprocess : str | dict | list | None
        Postprocessing config. See below for details.

    outfile : str
        Output filename to save final output. Expected to be gzipped.

    log : str
        Log filename that will accumulate all logs

    Notes
    -----

    This function:

        - downloads the URL[s] to tempfile[s]
        - resolves the name of the postprocessing function(s) if provided and
          imports it
        - calls the imported postprocessing function using the tempfile[s] and
          outfile plus any additional specified arguments.

    The postprocessing function must have one of the following signatures,
    where `infiles` contains the list of temporary files downloaded from the
    URL or URLs specified, and `outfile` is a gzipped file expected to be
    created by the function::

        def func(infiles, outfile):
            pass

    or::

        def func(infiles, outfile, *args):
            pass

    or::

        def func(infiles, outfile, *args, **kwargs):
            pass


    The function is specified as a string that resolves to an importable
    function, e.g., `postprocess: lib.postprocess.dm6.fix` will call a function
    called `fix` in the file `lib/postprocess/dm6.py`.

    If the contents of `postprocess:` is a dict, it must have at least the key
    `function`, and optionally `args` and/or `kwargs` keys. The `function` key
    indicates the importable path to the function.  `args` can be a string
    or list of arguments that will be provided as additional args to a function
    with the second kind of signature above.  If `kwargs` is provided, it is
    a dict that is passed to the function with the third kind of signature
    above. For example::

        postprocess:
            function: lib.postprocess.dm6.fix
            args:
                - True
                - 3

    or::

        postprocess:
            function: lib.postprocess.dm6.fix
            args:
                - True
                - 3
            kwargs:
                skip: exon

    """

    def default_postprocess(origfn, newfn):
        shell("mv {origfn} {newfn}")

    if not isinstance(postprocess, list):
        postprocess = [postprocess]

    # Will contain tuples of (func, args, kwargs, tmp_outfile)
    funcs = []

    # It is possible to chain multiple postprocessing functions together by
    # providing them as a list.
    #
    #   postprocess = [
    #
    #     "lib.func1",
    #
    #     {
    #       "function": "lib.func2",
    #       "args": (True, True),
    #     },
    #
    #     {
    #       "function": "lib.func3",
    #       "args": (1, 2),
    #       "kwargs": {"gzipped": True),
    #     },
    #
    #   ]
    #
    for i, postprocess_i in enumerate(postprocess):

        if postprocess_i is None:
            func = default_postprocess
            args = ()
            kwargs = {}
            name = None

        # postprocess can have a single string value indicating the function or
        # it can be a dict with keys "function" and optionally "args". The value of
        # "args" can be a string or a list.
        else:
            if isinstance(postprocess_i, dict):
                name = postprocess_i.get("function", postprocess)
                args = postprocess_i.get("args", ())
                kwargs = postprocess_i.get("kwargs", {})
                if isinstance(args, str):
                    args = (args,)
            elif isinstance(postprocess_i, str):
                name = postprocess_i
                args = ()
                kwargs = {}

            else:
                raise ValueError(
                    f"Unhandled type of postprocessing configuration: {postprocess_i}"
                )

            # In the special case where there is kwarg beginning and ending
            # with "__", this can be a dotted function name so it will be
            # resolved here as well and passed along to the postprocessing
            # function.
            #
            # This makes it possible to do things like add ERCC annotations on
            # the end of other annotations that themselves need to be
            # post-processed.
            for kw in kwargs:
                if kw.startswith("__") and kw.endswith("__"):
                    kwargs[kw] = resolve_name(kwargs[kw])

            # import the function
            func = resolve_name(name)

        tmp_outfile = f"{outfile}.{i}.{name}.tmp"
        funcs.append([func, args, kwargs, tmp_outfile])

    # The last func's outfile should be the final outfile
    funcs[-1][-1] = outfile

    # as described in the docstring above, functions are to assume a list of
    # urls
    if isinstance(urls, str):
        urls = [urls]

    # Download into reasonably-named temp filenames
    downloaded_tmpfiles = [f"{outfile}.{i}.tmp" for i in range(len(urls))]

    # For the first postprocess, its input will be all the downloaded files.
    postprocess_input = downloaded_tmpfiles
    try:
        # Copy (if local URI) or download into the specified temp files
        for url, tmpfile in zip(urls, downloaded_tmpfiles):
            if url.startswith("file:"):
                url = url.replace("file://", "")
                shell("cp {url} {tmpfile} 2> {log}")
            else:
                shell("wget {url} -O- > {tmpfile} 2> {log}")

        for func, args, kwargs, tmp_outfile in funcs:
            func(
                # all downloaded files (if the first postprocess), or the
                # output of the last postprocess
                postprocess_input,
                # the temp output for just this postprocess
                tmp_outfile,
                *args,
                **kwargs,
            )

            # We want the next postprocess to use the output of what we just
            # ran; as documented above the input files are expected to be in
            # a list.
            postprocess_input = [tmp_outfile]

    except Exception as e:
        raise e
    finally:
        to_delete = downloaded_tmpfiles

        # all but the last postprocess func output (the last one is the final
        # output that we want to keep!)
        to_delete += [i[-1] for i in funcs[:-1]]

        for i in to_delete:
            if os.path.exists(i):
                shell("rm {i}")
    if not _is_gzipped(outfile):
        raise ValueError(f"{outfile} does not appear to be gzipped.")


def get_sampletable(config):
    """
    Return samples and pandas.DataFrame of parsed sampletable.

    Returns the sample IDs and the parsed sampletable from the file specified
    in the config.

    The sample IDs are assumed to be the first column of the sampletable.

    Parameters
    ----------
    config : dict
    """
    sampletable = pandas.read_csv(config["sampletable"], comment="#", sep="\t")
    samples = sampletable.iloc[:, 0]
    return samples, sampletable


def get_techreps(sampletable, label):
    """
    Return all sample IDs for which the "label" column is `label`.
    """
    # since we're not requiring a name but we want to use `loc`
    first_col = sampletable.columns[0]
    result = list(sampletable.loc[sampletable["label"] == label, first_col])

    # If we're using a ChIP-seq-like sampletable we can provide a more
    # informative error message.

    is_chipseq = "antibody" in sampletable.columns
    if is_chipseq:
        err = """
        No technical replicates found for label '{}'. Check the ChIP-seq config
        file to ensure the peak-calling section only specifies values from the
        sampletable's "label" column.""".format(
            label
        )
    else:
        err = "No technical replicates found for label '{}'.".format(label)

    if len(result) == 0:
        raise ValueError(err)

    return result


def deprecation_handler(config):
    """
    Checks the config to see if anything has been deprecated.

    Also makes any fixes that can be done automatically.
    """
    if "assembly" in config:
        config["organism"] = config["assembly"]
        warnings.warn(
            "'assembly' should be replaced with 'organism' in config files. "
            "As a temporary measure, a new 'organism' key has been added with "
            "the value of 'assembly'",
            DeprecationWarning,
        )

    for org, block1 in config.get("references", {}).items():
        for tag, block2 in block1.items():
            gtf_conversions = block2.get("gtf", {}).get("conversions", [])
            for c in gtf_conversions:
                if isinstance(c, dict) and "annotation_hub" in c:
                    warnings.warn(
                        "You may want to try the 'mappings' conversion rather "
                        "than 'annotation_hub' since it works directly off "
                        "the GTF file rather than assuming concordance between "
                        "GTF and AnnoationHub instances",
                        DeprecationWarning,
                    )

    return config


def check_url(url, verbose=False):
    """
    Try to open -- and then immediately close -- a URL.

    Any exceptions can be handled upstream.

    """

    # Some notes here:
    #
    #  - A pure python implementation isn't great because urlopen seems to
    #    cache or hold sessions open or something. EBI servers reject responses
    #    because too many clients are connected. This doesn't happen using curl.
    #
    #  - Using the requests module doesn't help, because urls can be ftp:// and
    #    requests doesn't support that.
    #
    #  - Similarly, using asyncio and aiohttp works great for https, but not
    #    ftp (I couldn't get aioftp to work properly).
    #
    #  - Not all servers support --head. An example of this is
    #    https://www-s.nist.gov/srmors/certificates/documents/SRM2374_Sequence_v1.FASTA.
    #
    #  - Piping curl to head using the -c arg to use bytes seems to work.
    #    However, we need to set pipefail (otherwise because head exits 0 the
    #    whole thing exits 0). And in that case, we expect curl to exit every
    #    time with exit code 23, which is "failed to write output", because of
    #    the broken pipe. This is handled below.
    #
    if verbose:
        print(f"Checking {url}")

    # Notes on curl args:
    #
    #  --max-time to allow the server some seconds to respond
    #  --retry to allow multiple tries if transient errors (4xx for FTP, 5xx for HTTP) are found
    #  --silent to not print anything
    #  --fail to return non-zero exit codes for 404 (default is exit 0 on hitting 404)
    #
    # Need to run through bash explicitly to get the pipefail option, which in
    # turn means running with shell=True
    proc = subprocess.run(
        f'/bin/bash -o pipefail -c "curl --retry 3 --max-time 10 --silent --fail {url} | head -c 10 > /dev/null"',
        shell=True,
    )
    return proc


def check_urls(config, verbose=False):
    """
    Given a config filename or existing object, extract the URLs and check
    them.

    Parameters
    ----------

    config : str or dict
        Config object to inspect

    verbose : bool
        Print which URL is being checked

    wait : int
        Number of seconds to wait in between checking URLs, to avoid
        too-many-connection issues
    """
    failures = []
    urls = list(set(utils.flatten(pluck(config, "url"))))
    for url in urls:
        if url.startswith("file://"):
            continue

        res = check_url(url, verbose=verbose)

        # we expect exit code 23 because we're triggering SIGPIPE with the
        # "|head -c" above.
        if res.returncode and res.returncode != 23:
            failures.append(
                f"FAIL with exit code {res.returncode}. Command was: {res.args}"
            )
    if failures:
        output = "\n   ".join(failures)
        raise ValueError(
            f"Found problematic URLs. See https://ec.haxx.se/usingcurl/usingcurl-returns for explanation of exit codes.\n   {output}"
        )


def check_all_urls_found(verbose=True):
    """
    Recursively loads all references that can be included and checks them.
    Reports out if there are any failures.
    """
    check_urls(
        {
            "include_references": [
                "include/reference_configs",
                "test/test_configs",
                "workflows/rnaseq/config",
                "workflows/chipseq/config",
                "workflows/references/config",
            ]
        },
        verbose=verbose,
    )


def gff2gtf(gff, gtf):
    """
    Converts a gff file to a gtf format using the gffread function from Cufflinks
    """
    if _is_gzipped(gff[0]):
        shell("gzip -d -S .gz.0.tmp {gff} -c | gffread - -T -o- | gzip -c > {gtf}")
    else:
        shell("gffread {gff} -T -o- | gzip -c > {gtf}")


def wrapper_for(path):
    return 'file:' + os.path.join('../..','wrappers', 'wrappers', path)

def detect_sra(sampletable):
    return 'Run' in self.sampletable.columns and any(self.sampletable['Run'].str.startswith('SRR'))

# vim: ft=python
