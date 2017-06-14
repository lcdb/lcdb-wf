import os
import tempfile
import yaml
import pandas
from Bio import SeqIO
import gzip
from lcdblib.utils.imports import resolve_name
from lcdblib.snakemake import aligners
from snakemake.shell import shell


def gzipped(tmpfiles, outfile):
    """
    Cat-and-gzip input files into a single output file.
    """
    with gzip.open(outfile, 'wt') as fout:
        for f in tmpfiles:
            for line in open(f):
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
    Given input FASTAs, create a new one containing only records whose
    description matches `pattern`.

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


def download_and_postprocess(outfile, config, assembly, tag, type_):
    """
    Given an output file, figure out what to do based on the config.

    Parameters
    ----------
    outfile : str

    config : dict

    assembly : str
        Which assembly to use. Must be a key in the "references" section of the
        config.

    tag : str
        Which tag for the assembly to use. Must be a tag for the assembly in the config

    type_ : str
        A supported references type (gtf, fasta) to use.

    Notes
    -----

    This function:

     - uses `assembly`, `tag`, `type_` as a key into the config dict to figure
       out:
         - what postprocessing function (if any) was specified along with
           its optional args
         - the URL[s] to download
     - resolves the name of that function and imports it
     - downloads the URL[s] to tempfile[s]
     - calls the imported function using the tempfile[s] and outfile plus
       any additional specified arguments.


    The function must have one of the following two signatures::

        def func(infiles, outfile):
            pass

    or::

        def func(infiles, outfile, *args):
            pass


    `infiles` contains the list of temporary files downloaded from the URL or
    URLs specified.

    `outfile` is a gzipped file expected to be created by the function.

    The function is specified as a string that resolves to an importable
    function, e.g., `lib.postprocess.dm6.fix` is a function called `fix` in the
    file `lib/postprocess/dm6.py`.

    If specified in the config as a dict, it must have `function` and `args`
    keys. The `function` key indicates the importable path to the function, and
    `args` can be a string or list of arguments that will be provided as
    additional args to a function with the second kind of signature above.
    """

    def default_postprocess(origfn, newfn):
        """
        If no other postprocess function is defined, then simply move the original
        to the new.
        """
        shell("mv {origfn} {newfn}")

    block = config['references'][assembly][tag][type_]

    # postprocess can be missing, in which case we use the default above
    post_process = block.get('postprocess', None)
    if post_process is None:
        func = default_postprocess
        args = ()

    # postprocess can have a single string value (indicating the function) or
    # it can be a dict with keys "function" and optionally "args". The value of
    # "args" can be a string or a list.
    else:
        if isinstance(post_process, dict):
            name = post_process.get('function', post_process)
            args = post_process.get('args', ())
            if isinstance(args, str):
                args = (args,)
        elif isinstance(post_process, str):
            name = post_process
            args = ()

        # import the function
        func = resolve_name(name)

    # functions assume a list of urls
    urls = block['url']
    if isinstance(urls, str):
        urls = [urls]

    tmpfiles = ['{0}.{1}.tmp'.format(outfile, i) for i in range(len(urls))]
    try:
        for url, tmpfile in zip(urls, tmpfiles):
            shell("wget {url} -O- > {tmpfile} 2> {outfile}.log")

        func(tmpfiles, outfile, *args)
    except Exception as e:
        raise e
    finally:
        for i in tmpfiles:
            if os.path.exists(i):
                shell('rm {i}')


def references_dict(config):
    """
    Reformats the config file's reference section into a more practical form.

    Files can be referenced as `d[assembly][tag][type]`.

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
    ...           - kallisto
    ... '''))
    >>> fout.close()

    To this format:

    >>> d = config_to_dict('tmp')
    >>> assert d == (
    ... {
    ...   'dm6': {
    ...      'r6-11': {
    ...          'fasta': '/data/dm6/r6-11/fasta/dm6_r6-11.fasta',
    ...          'refflat': '/data/dm6/r6-11/gtf/dm6_r6-11.refflat',
    ...          'gtf': '/data/dm6/r6-11/gtf/dm6_r6-11.gtf',
    ...          'chromsizes': '/data/dm6/r6-11/fasta/dm6_r6-11.chromsizes',
    ...          'bowtie2': '/data/dm6/r6-11/bowtie2/dm6_r6-11.1.bt2',
    ...          'hisat2': '/data/dm6/r6-11/hisat2/dm6_r6-11.1.ht2',
    ...          },
    ...      'r6-11_transcriptome': {
    ...          'fasta': '/data/dm6/r6-11_transcriptome/fasta/dm6_r6-11_transcriptome.fasta',
    ...          'chromsizes': '/data/dm6/r6-11_transcriptome/fasta/dm6_r6-11_transcriptome.chromsizes',
    ...          'kallisto': '/data/dm6/r6-11_transcriptome/kallisto/dm6_r6-11_transcriptome.idx',
    ...          },
    ...     },
    ... }), d
    >>> os.unlink('tmp')

    """
    if isinstance(config, str):
        config = yaml.load(open(config))

    references_dir = get_references_dir(config)

    # Map "indexes" value to a pattern specific to each index.
    index_extensions = {
        'bowtie2': aligners.bowtie2_index_from_prefix('')[0],
        'hisat2': aligners.hisat2_index_from_prefix('')[0],
        'kallisto': '.idx',
        'salmon': '/hash.bin'
    }

    conversion_extensions = {
        'intergenic': '.intergenic.gtf',
        'refflat': '.refflat',
        'gffutils': '.gtf.db',
        'genelist': '.genelist',
        'annotation_hub': '.{keytype}.csv'
    }

    d = {}
    conversion_kwargs = {}
    for assembly in config['references'].keys():
        d[assembly] = {}
        for tag in config['references'][assembly].keys():
            e = {}
            for type_, block in config['references'][assembly][tag].items():
                e[type_] = (
                    '{references_dir}/'
                    '{assembly}/'
                    '{tag}/'
                    '{type_}/'
                    '{assembly}_{tag}.{type_}'.format(**locals())
                )

                # Add conversions if specified.
                if type_ == 'gtf':
                    conversions = block.get('conversions', [])
                    for conversion in conversions:
                        kwargs = {}
                        if isinstance(conversion, dict):
                            # we assume that there is only one key, and that
                            # key is the actual name of the conversion; the
                            # corresponding value will be kwargs
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
                            '{assembly}/'
                            '{tag}/'
                            '{type_}/'
                            '{assembly}_{tag}{ext}'.format(**locals())
                        )
                        e[conversion] = output

                        conversion_kwargs[output] = kwargs

                if type_ == 'fasta':
                    # Add indexes if specified
                    indexes = block.get('indexes', [])
                    for index in indexes:
                        ext = index_extensions[index]

                        e[index] = (
                            '{references_dir}/{assembly}/{tag}/{index}/{assembly}_{tag}{ext}'
                            .format(**locals())
                        )

                    e['chromsizes'] = (
                        '{references_dir}/'
                        '{assembly}/'
                        '{tag}/'
                        '{type_}/'
                        '{assembly}_{tag}.chromsizes'.format(**locals())
                    )
                d[assembly][tag] = e
    return d, conversion_kwargs


def tempdir_for_biowulf():
    """
    Get an appropriate tempdir.

    The NIH biowulf cluster allows nodes to have their own /lscratch dirs as
    local temp storage. However the particular dir depends on the slurm job ID,
    which is not known in advance. This function sets the shell.prefix to use
    such a tempdir if it's available; otherwise it leaves TMPDIR unchanged.
    This makes it suitable for running locally or on other clusters, however if
    you need different behavior then a different function will need to be
    written.
    """
    tmpdir = tempfile.gettempdir()
    jobid = os.getenv('SLURM_JOBID')
    if jobid:
        tmpdir = os.path.join('/lscratch', jobid)
    return tmpdir


def get_references_dir(config):
    """
    Returns the references dir, preferring the value of an existing environment
    variable `REFERENCES_DIR` over the config entry "references_dir". Raise an
    error if either can't be found.

    Parameters
    ----------
    config : dict
    """
    references_dir = os.environ.get(
        'REFERENCES_DIR', config.get('references_dir', None))
    if references_dir is None:
        raise ValueError('No references dir specified')
    config['references_dir'] = references_dir
    return references_dir


def get_sampletable(config):
    """
    Returns the sample IDs and the parsed sampletable.

    The sample IDs are assumed to be the first column of the sampletable.

    Parameters
    ----------
    config : dict
    """
    sampletable = pandas.read_table(config['sampletable'])
    samples = sampletable.iloc[:, 0]
    return samples, sampletable
