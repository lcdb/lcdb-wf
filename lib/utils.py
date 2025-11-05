import binascii
import csv
import gzip
import os
import subprocess
import warnings
from collections.abc import Iterable

import gffutils
import pandas as pd
from Bio import SeqIO
from snakemake.shell import shell


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
            module_ = __import__(".".join(parts_copy))
            break
        except ImportError:
            del parts_copy[-1]
            if not parts_copy:
                raise
    parts = parts[1:]
    obj = module_
    for part in parts:
        obj = getattr(obj, part)
    return obj


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


def is_gzipped(fn):
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
    if is_gzipped(tmp):
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


def is_paired_end(sampletable, sample):
    """
    Inspects the sampletable to see if the sample is paired-end or not

    Parameters
    ----------
    sampletable : pandas.DataFrame
        If SRA sampletable, contains a "layout" or "LibraryLayout" column (but
        not both). If the lowercase value is "pe" or "paired", consider the
        sample paired-end. Otherwise consider single-end.

        Otherwise, if there's an "orig_filename_R2" column consider it
        paired-end, otherwise single-end.

    sample : str
        Assumed to be found in the first column of `sampletable`
    """
    # We can't fall back to detecting PE based on two fastq files provided for
    # each sample when it's an SRA sampletable (which only has SRR accessions).
    #
    # So detect first detect if SRA sampletable based on presence of "Run"
    # column and all values of that column starting with "SRR", and then raise
    # an error if the Layout or LibraryLayout column does not exist.

    sra_layout_columns = ["layout", "librarylayout"]
    sampletable_columns = [i.lower() for i in sampletable.columns]
    if "run" in sampletable_columns:
        if all(sampletable["Run"].str.startswith("SRR")):
            if len(set(sra_layout_columns).intersection(sampletable_columns)) == 0:
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
            report_ = f"SE samples: {s}"
        else:
            report_ = f"PE samples: {p}"
        raise ValueError(f"Only a single layout (SE or PE) is supported. {report_}")


def filter_rrna_fastas(tmpfiles, outfile, pattern):
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
    if pattern is None:
        raise ValueError("Pattern cannot be None")

    def gen():
        for tmp in tmpfiles:
            handle = gzip.open(tmp, "rt")
            parser = SeqIO.parse(handle, "fasta")
            for rec in parser:
                if pattern not in rec.description:
                    continue
                rec.seq = rec.seq.back_transcribe()
                # rec.description = rec.name
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
    def _default(origfn, newfn):
        shell("mv {origfn} {newfn}")

    for i, postprocess_i in enumerate(postprocess):

        if postprocess_i is None:
            func = _default
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
    if not is_gzipped(outfile):
        raise ValueError(f"{outfile} does not appear to be gzipped.")


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
    urls = list(set(flatten(pluck(config, "url"))))
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


def gff2gtf(gff, gtf):
    """
    Converts a gff file to a gtf format using the gffread function from Cufflinks
    """
    if is_gzipped(gff[0]):
        shell("gzip -d -S .gz.0.tmp {gff} -c | gffread - -T -o- | gzip -c > {gtf}")
    else:
        shell("gffread {gff} -T -o- | gzip -c > {gtf}")


def detect_sra(sampletable):
    return "Run" in sampletable.columns and any(
        sampletable["Run"].str.startswith("SRR")
    )


def mappings_tsv(
    gtf,
    tsv,
    exclude_featuretypes=None,
    include_featuretypes=None,
    include_attributes=None,
    exclude_attributes=None,
):
    """
    Create a TSV file of attributes found in a GTF file.

    Parameters
    ----------

    gtf, tsv : str
        Input and output filenames respectively

    exclude_featuretypes, include_featuretypes : list
        Mutually exclusive; use these to restrict the features considered.
        E.g., we likely don't need entries for start_codon if those are in the
        GTF.

    include_attributes, exclude_attributes : list
        Mutually exclusive. Restrict the attributes reported in the TSV. Should at least have
        a column for gene ID and transcript ID in order for downstream RNA-seq
        work.
    """

    if exclude_featuretypes and include_featuretypes:
        raise ValueError(
            "Both include_featuretypes and exclude_featuretypes were specified."
        )
    if exclude_attributes and include_attributes:
        raise ValueError(
            "Both include_attributes and exclude_attributes were specified."
        )

    res = []
    keys = set(["__featuretype__"])
    seen = set()
    for f in gffutils.DataIterator(gtf):
        ft = f.featuretype
        if exclude_featuretypes and ft in exclude_featuretypes:
            continue
        if include_featuretypes and ft not in include_featuretypes:
            continue

        d = dict(f.attributes)

        if include_featuretypes:
            d = {k: v for k, v in d.items() if k in include_featuretypes}
        if exclude_featuretypes:
            d = {k: v for k, v in d.items() if k not in exclude_featuretypes}

        keys.update(d.keys())
        d["__featuretype__"] = ft

        # Exclude duplicates (rather than sorting and uniq-ing the file later)
        h = hash(str(d))
        if h in seen:
            continue
        seen.update([h])

        res.append(d)

    def unlist_dict(d):
        for k, v in d.items():
            if isinstance(v, list):
                d[k] = "|".join(v)
        return d

    if include_attributes:
        sorted_keys = sorted(include_attributes)
    else:
        sorted_keys = sorted(keys)
    with open(tsv, "w") as fout:
        writer = csv.DictWriter(
            fout, fieldnames=sorted_keys, restval="", delimiter="\t"
        )
        writer.writeheader()
        for row in res:
            writer.writerow(unlist_dict(row))


def preflight(config, sampletable):
    """
    Performs verifications on config and sampletable files

    Parameters
    ----------
    config: yaml config object
    """

    if len(sampletable) != len(sampletable.iloc[:, 0].unique()):
        raise ConfigurationError("Samplenames non unique, check the sampletable")

    # For non-SRA sampletables
    if "orig_filename" in sampletable.columns:
        fns = sampletable["orig_filename"]
        if "orig_filename_R2" in sampletable.columns:
            fns = pd.concat([fns, sampletable["orig_filename_R2"]])
        if len(fns.unique()) < len(fns):
            raise ValueError("Fastq filenames non unique, check the sampletable\n")

    if "genome" not in config:
        raise ConfigurationError("Config is missing 'genome' key")
    if "url" not in config["genome"]:
        raise ConfigurationError("Config is missing 'url' key for 'genome'")


def rnaseq_preflight(config, sampletable):
    preflight(config, sampletable)
    if "annotation" not in config:
        raise ConfigurationError("Config is missing 'annotation' key")
    if "url" not in config["annotation"]:
        raise ConfigurationError("Config is missing 'url' key for 'annotation'")
    if "stranded" not in config:
        raise ConfigurationError("Config is missing 'stranded' key")
    if "organism" not in config:
        raise ConfigurationError("Config is missing 'organism' key")


def chipseq_preflight(config, sampletable):
    preflight(config, sampletable)
    if "peaks" not in config:
        config["peaks"] = []


def read_sampletable(config):
    """
    Given a config object, return the sampletable with the first column used as the index.

    Autodetect tsv/csv.
    """
    sampletable_fn = config.get("sampletable", "config/sampletable.tsv")
    if sampletable_fn.endswith(".tsv"):
        sep = "\t"
    elif sampletable_fn.endswith(".csv"):
        sep = ","
    else:
        raise ConfigurationError(
            f"Sampletable should end in .csv or .tsv to indicate format, got {sampletable_fn}"
        )
    sampletable = pd.read_table(sampletable_fn, sep=sep, comment="#")
    sampletable = sampletable.set_index(sampletable.columns[0], drop=False)
    return sampletable


def prepare_chipseq_sampletable(config):
    """
    Given a config, return the validated and prepared ChIP-seq table.
    """
    sampletable = read_sampletable(config)
    sampletable["label"] = sampletable["label"].fillna(sampletable.iloc[:, 0])
    chipseq_preflight(config, sampletable)
    return sampletable


def prepare_rnaseq_sampletable(config):
    """
    Given a config, return the validated and prepared RNA-seq table.
    """
    sampletable = read_sampletable(config)
    rnaseq_preflight(config, sampletable)
    return sampletable


# vim: ft=python
