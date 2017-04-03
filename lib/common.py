import os
import yaml
from Bio import SeqIO
import gzip
from lcdblib.utils.imports import resolve_name
from lcdblib.snakemake import aligners, helpers
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

    This function:

     - uses assembly, tag, type_ as a key into the config dict to figure out:
         - what postprocessing function (if any) was specified along with
           its optional args
         - the URL[s] to download
     - resolves the name of that function and imports it
     - downloads the URL[s] to tempfile[s]
     - calls the imported function using the tempfile[s] and outfile plus
       any additional specified arguments.

    If defined, the function must assume a list of input gzipped files and must
    create the gzipped output file (whose name is given as its 2nd input arg).
    """
    references_dir = get_references_dir(config)

    def default_postprocess(origfn, newfn):
        """
        If no other postprocess function is defined, then simply move the original
        to the new.
        """
        shell("mv {origfn} {newfn}")

    base = os.path.relpath(outfile, references_dir)
    basename = os.path.basename(outfile)

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


def get_references_dir(config):
    references_dir = os.environ.get('REFERENCES_DIR', config.get('references_dir', None))
    if references_dir is None:
        raise ValueError('References dir not set')
    return references_dir


def references_dict(config):
    """
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

    With this new dictionary, other parts of the config or snakemake rules can
    access the files by d[assembly][tag][type].
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

                if type_== 'fasta':
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
