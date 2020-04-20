import glob
import os
import warnings
import yaml
import pandas
from Bio import SeqIO
import gzip
import binascii
from lib.imports import resolve_name
from lib import aligners
from snakemake.shell import shell
from snakemake.io import expand

# List of possible keys in config that are to be interpreted as paths
PATH_KEYS = [
    'references_dir',
    'sampletable',
    'sample_dir',
    'aggregation_dir',
    'merged_dir',
    'peaks_dir',
    'hub_config',
]


def _is_gzipped(fn):
    """
    Filename-independent method of checking if a file is gzipped or not. Uses
    the magic number.

    xref https://stackoverflow.com/a/47080739
    """
    with open(fn, 'rb') as f:
        return binascii.hexlify(f.read(2)) == b'1f8b'


def openfile(tmp, mode):
    """
    Returns an open file handle; auto-detects gzipped files.
    """
    if _is_gzipped(tmp):
        return gzip.open(tmp, mode)
    else:
        return open(tmp, mode)


def resolve_config(config, workdir=None):
    """
    Finds the config file.

    Parameters
    ----------
    config : str, dict
        If str, assume it's a YAML file and parse it; otherwise pass through

    workdir : str
        Optional location to specify relative location of all paths in `config`
    """
    if isinstance(config, str):
        config = yaml.load(open(config), Loader=yaml.FullLoader)

    def rel(pth):
        if workdir is None or os.path.isabs(pth):
            return pth
        return os.path.join(workdir, pth)
    for key in PATH_KEYS:
        if key in config:
            config[key] = rel(config[key])
    return config


def gzipped(tmpfiles, outfile):
    """
    Cat-and-gzip a list of uncompressed files into a compressed output file.
    """
    with gzip.open(outfile, 'wt') as fout:
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
    shell('cat {tmpfiles} > {outfile}')


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
            handle = gzip.open(tmp, 'rt')
            parser = SeqIO.parse(handle, 'fasta')
            for rec in parser:
                if pattern not in rec.description:
                    continue
                rec.seq = rec.seq.back_transcribe()
                rec.description = rec.name
                yield rec

    with gzip.open(outfile, 'wt') as fout:
        SeqIO.write(gen(), fout, 'fasta')


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
    lookup = {i: i + '.fa' for i in tmpfiles}
    for i in tmpfiles:
        fn = lookup[i]
        shell('twoBitToFa {i} {fn}')

    # Make sure we retain the order of the originally-provided files from the
    # config when concatenating.
    fastas = [lookup[i] for i in tmpfiles]
    shell('cat {fastas} | gzip -c > {outfile}')
    shell('rm {fastas}')


def download_and_postprocess(outfile, config, organism, tag, type_):
    """
    Given an output file, figure out what to do based on the config.

    See notes below for details.

    Parameters
    ----------
    outfile : str

    config : dict

    organism : str
        Which organism to use. Must be a key in the "references" section of the
        config.

    tag : str
        Which tag for the organism to use. Must be a tag for the organism in
        the config

    type_ : str
        A supported references type (gtf, fasta) to use.

    Notes
    -----

    This function:

        - uses `organism`, `tag`, `type_` as a key into the config dict to
          figure out:

            - what postprocessing function (if any) was specified along with
              its optional args
            - the URL[s] to download

        - resolves the name of the postprocessing function (if provided) and
          imports it
        - downloads the URL[s] to tempfile[s]
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
        """
        If no other postprocess function is defined, then simply move the
        original to the new.
        """
        shell("mv {origfn} {newfn}")

    block = config['references'][organism][tag][type_]

    # postprocess can be missing, in which case we use the default above
    post_process = block.get('postprocess', None)

    if not isinstance(post_process, list):
        post_process = [post_process]

    funcs = []
    func_tmpfiles = []
    for i, post_process_block in enumerate(post_process):
        if post_process_block is None:
            func = default_postprocess
            args = ()
            kwargs = {}
            name = None

        # postprocess can have a single string value (indicating the function) or
        # it can be a dict with keys "function" and optionally "args". The value of
        # "args" can be a string or a list.
        else:
            if isinstance(post_process_block, dict):
                name = post_process_block.get('function', post_process)
                args = post_process_block.get('args', ())
                kwargs = post_process_block.get('kwargs', {})
                if isinstance(args, str):
                    args = (args,)
            elif isinstance(post_process_block, str):
                name = post_process_block
                args = ()
                kwargs = {}

            # import the function
            func = resolve_name(name)

        tmp_outfile = f'{outfile}.{i}.{name}.tmp'
        func_tmpfiles.append(tmp_outfile)
        funcs.append([func, args, kwargs, tmp_outfile])

    # The last func's outfile should be the final outfile
    funcs[-1][-1] = outfile

    # as described in the docstring above, functions are to assume a list of
    # urls
    urls = block['url']
    if isinstance(urls, str):
        urls = [urls]

    # Download tempfiles into reasonably-named filenames
    tmpfiles = ['{0}.{1}.tmp'.format(outfile, i) for i in range(len(urls))]
    tmpinputfiles = tmpfiles
    try:
        for url, tmpfile in zip(urls, tmpfiles):
            if url.startswith('file:'):
                url = url.replace('file://', '')
                shell('cp {url} {tmpfile} 2> {outfile}.log')
            else:
                shell("wget {url} -O- > {tmpfile} 2> {outfile}.log")

        for func, args, kwargs, outfile in funcs:
            func(tmpinputfiles, outfile, *args, **kwargs)
            tmpinputfiles = [outfile]

    except Exception as e:
        raise e
    finally:
        for i in tmpfiles + func_tmpfiles:
            if os.path.exists(i):
                shell('rm {i}')


def references_dict(config):
    """
    Reformats the config file's reference section into a more practical form.

    Files can be referenced as `d[organism][tag][type]`.

    Parameters
    ----------
    config : dict

    Notes
    -----

    The config file is designed to be easy to edit and use from the user's
    standpoint. But it's not so great for practical usage. Here we convert the
    config file which has the format::

    >>> from textwrap import dedent
    >>> fout = open('tmp', 'w')
    >>> _ = fout.write(dedent('''
    ... references_dir: "/data"
    ... references:
    ...   dm6:
    ...     r6-11:
    ...       metadata:
    ...         reference_genome_build: 'dm6'
    ...         reference_effective_genome_count: 1.2e7
    ...         reference_effective_genome_proportion: 0.97
    ...       fasta:
    ...         url: ""
    ...         indexes:
    ...           - bowtie2
    ...           - hisat2
    ...       gtf:
    ...         url: ""
    ...         conversions:
    ...           - refflat
    ...     r6-11_transcriptome:
    ...       fasta:
    ...         url: ""
    ...         indexes:
    ...           - salmon
    ... '''))
    >>> fout.close()

    To this format:

    >>> d, conversion_kwargs = references_dict('tmp')
    >>> assert d == (
    ... {
    ...   'dm6': {
    ...      'r6-11': {
    ...          'fasta': '/data/dm6/r6-11/fasta/dm6_r6-11.fasta',
    ...          'refflat': '/data/dm6/r6-11/gtf/dm6_r6-11.refflat',
    ...          'gtf': '/data/dm6/r6-11/gtf/dm6_r6-11.gtf',
    ...          'chromsizes': '/data/dm6/r6-11/fasta/dm6_r6-11.chromsizes',
    ...          'bowtie2': '/data/dm6/r6-11/bowtie2/dm6_r6-11.1.bt2',
    ...          'bowtie2_fasta': '/data/dm6/r6-11/bowtie2/dm6_r6-11.fasta',
    ...          'hisat2': '/data/dm6/r6-11/hisat2/dm6_r6-11.1.ht2',
    ...          'hisat2_fasta': '/data/dm6/r6-11/hisat2/dm6_r6-11.fasta',
    ...          },
    ...      'r6-11_transcriptome': {
    ...          'fasta': '/data/dm6/r6-11_transcriptome/fasta/dm6_r6-11_transcriptome.fasta',
    ...          'chromsizes': '/data/dm6/r6-11_transcriptome/fasta/dm6_r6-11_transcriptome.chromsizes',
    ...          'salmon': '/data/dm6/r6-11_transcriptome/salmon/dm6_r6-11_transcriptome/hash.bin',
    ...          'salmon_fasta': '/data/dm6/r6-11_transcriptome/salmon/dm6_r6-11_transcriptome.fasta',
    ...          },
    ...     },
    ... }), d
    >>> assert conversion_kwargs == {'/data/dm6/r6-11/gtf/dm6_r6-11.refflat': {}}
    >>> os.unlink('tmp')

    """
    if isinstance(config, str):
        config = yaml.load(open(config), Loader=yaml.FullLoader)

    references_dir = get_references_dir(config)

    # Map "indexes" value to a pattern specific to each index.
    index_extensions = {
        'bowtie2': aligners.bowtie2_index_from_prefix('')[0],
        'hisat2': aligners.hisat2_index_from_prefix('')[0],
        'star': '/Genome',
        'salmon': '/hash.bin',
    }

    conversion_extensions = {

        'intergenic': '.intergenic.gtf',
        'refflat': '.refflat',
        'gffutils': '.gtf.db',
        'bed12': '.bed12',
        'genelist': '.genelist',
        'annotation_hub': '.{keytype}.csv',
        'mappings': '.mapping.tsv.gz',
    }

    d = {}
    conversion_kwargs = {}

    merged_references = config['references']

    for organism in merged_references.keys():
        d[organism] = {}
        for tag in merged_references[organism].keys():
            e = {}
            for type_, block in merged_references[organism][tag].items():
                if type_ == 'metadata':
                    continue
                e[type_] = (
                    '{references_dir}/'
                    '{organism}/'
                    '{tag}/'
                    '{type_}/'
                    '{organism}_{tag}.{type_}'.format(**locals())
                )

                # Add conversions if specified.
                if type_ == 'gtf':
                    conversions = block.get('conversions', [])
                    for conversion in conversions:
                        kwargs = {}
                        if isinstance(conversion, dict):
                            # if conversion is specified as dict, we assume
                            # that there is only one key, and that key is the
                            # actual name of the conversion; the corresponding
                            # value will be kwargs. This is used e.g. for
                            # gffutils conversion which often need some
                            # tweaking of args depending on the gtf format.
                            assert len(list(conversion.keys())) == 1
                            kwargs = list(conversion.values())[0]
                            conversion = list(conversion.keys())[0]

                        # While the full set of columns for annotation hub are
                        # not known in advance, we can assume at least the
                        # keytype provided will be an output file. Fill that in
                        # here.
                        if conversion == 'annotation_hub':
                            keytype = kwargs['keytype']
                            ext = conversion_extensions[conversion].format(keytype=keytype)
                        else:
                            ext = conversion_extensions[conversion]
                        output = (
                            '{references_dir}/'
                            '{organism}/'
                            '{tag}/'
                            '{type_}/'
                            '{organism}_{tag}{ext}'.format(**locals())
                        )
                        e[conversion] = output

                        conversion_kwargs[output] = kwargs

                if type_ == 'fasta':
                    # Add indexes if specified
                    indexes = block.get('indexes', [])
                    for index in indexes:
                        ext = index_extensions[index]

                        e[index] = (
                            '{references_dir}/{organism}/{tag}/{index}/{organism}_{tag}{ext}'
                            .format(**locals())
                        )

                        # Each index will get the original fasta symlinked over
                        # to its directory
                        e[index + '_fasta'] = (
                            '{references_dir}/{organism}/{tag}/{index}/{organism}_{tag}.fasta'
                            .format(**locals())
                        )

                    e['chromsizes'] = (
                        '{references_dir}/'
                        '{organism}/'
                        '{tag}/'
                        '{type_}/'
                        '{organism}_{tag}.chromsizes'.format(**locals())
                    )
                d[organism][tag] = e
    return d, conversion_kwargs


def get_references_dir(config):
    """
    Identify the references directory based on config and env vars.

    Returns the references dir, preferring the value of an existing environment
    variable `REFERENCES_DIR` over the config entry "references_dir". Raise an
    error if either can't be found.

    Parameters
    ----------
    config : dict
    """
    config = resolve_config(config)
    references_dir = os.environ.get(
        'REFERENCES_DIR', config.get('references_dir', None))
    if references_dir is None:
        raise ValueError('No references dir specified')
    return references_dir


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
    config = resolve_config(config)
    sampletable = pandas.read_csv(config['sampletable'], comment="#", sep='\t')
    samples = sampletable.iloc[:, 0]
    return samples, sampletable


def get_techreps(sampletable, label):
    """
    Return all sample IDs for which the "label" column is `label`.
    """
    # since we're not requiring a name but we want to use `loc`
    first_col = sampletable.columns[0]
    result = list(sampletable.loc[sampletable['label'] == label, first_col])

    # If we're using a ChIP-seq-like sampletable we can provide a more
    # informative error message.

    is_chipseq = 'antibody' in sampletable.columns
    if is_chipseq:
        err = (
            "No technical replicates found for label '{}'. This looks to "
            "be a ChIP-seq experiment; check the peak-calling section of"
            "the config.".format(label)
        )
    else:
        err = "No technical replicates found for label '{}'.".format(label)

    if len(result) == 0:
        raise ValueError(err)

    return result


def load_config(config):
    """
    Loads the config.

    Resolves any included references directories/files and runs the deprecation
    handler.
    """
    if isinstance(config, str):
        config = yaml.load(open(config), Loader=yaml.FullLoader)

    # Here we populate a list of reference sections. Items later on the list
    # will have higher priority
    includes = config.get('include_references', [])
    reference_sections = []

    # First the directories. Directories that come earlier lose to those that
    # come later.
    for dirname in filter(os.path.isdir, includes):
        # Note we're looking recursively for .yaml and .yml, so very large
        # reference directories are possible
        for fn in glob.glob(os.path.join(dirname, '**/*.y?ml'),
                            recursive=True):
            refs = yaml.load(open(fn), Loader=yaml.FullLoader).get('references', None)
            if refs is None:
                raise ValueError("No 'references:' section in {0}".format(fn))
            reference_sections.append(refs)

    # Now the files
    for fn in filter(os.path.isfile, includes):
        refs = yaml.load(open(fn), Loader=yaml.FullLoader).get('references', None)
        if refs is None:
            raise ValueError("No 'references:' section in {0}".format(fn))
        reference_sections.append(refs)

    # The last thing we include is the references section as written in the
    # config, which wins over all.
    reference_sections.append(config.get('references', {}))

    merged_references = {}
    for ref in reference_sections:
        for organism in ref.keys():
            org_dict = merged_references.get(organism, {})
            for tag in ref[organism].keys():
                org_dict[tag] = ref[organism][tag]
            merged_references[organism] = org_dict
    config['references'] = merged_references

    # Run the deprecation handler on the final config
    config = deprecation_handler(config)

    return config


def deprecation_handler(config):
    """
    Checks the config to see if anything has been deprecated.

    Also makes any fixes that can be done automatically.
    """
    if 'assembly' in config:
        config['organism'] = config['assembly']
        warnings.warn(
            "'assembly' should be replaced with 'organism' in config files. "
            "As a temporary measure, a new 'organism' key has been added with "
            "the value of 'assembly'",
            DeprecationWarning)

    for org, block1 in config.get('references', {}).items():
        for tag, block2 in block1.items():
            gtf_conversions = block2.get('gtf', {}).get('conversions', [])
            for c in gtf_conversions:
                if isinstance(c, dict) and 'annotation_hub' in c:
                    warnings.warn(
                        "You may want to try the 'mappings' conversion rather "
                        "than 'annotation_hub' since it works directly off "
                        "the GTF file rather than assuming concordance between "
                        "GTF and AnnoationHub instances",
                        DeprecationWarning)

    return config


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
    row = sampletable.set_index(sampletable.columns[0]).loc[sample]
    if 'orig_filename_R2' in row:
        return True
    if 'layout' in row and 'LibraryLayout' in row:
        raise ValueError("Expecting column 'layout' or 'LibraryLayout', "
                         "not both")
    try:
        return row['layout'].lower() in ['pe', 'paired']
    except KeyError:
        pass
    try:
        return row['LibraryLayout'].lower() in ['pe', 'paired']
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
                'Need "{{sample}}" in pattern '
                '"{pattern}"'.format(pattern=pattern))
        n = [1]
        if is_paired_end(sampletable, wc.sample) and not r1_only:
            n = [1, 2]
        res = expand(pattern, sample=wc.sample, n=n)
        return res
    return func
