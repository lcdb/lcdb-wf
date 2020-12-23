"""
Helper functions for working with aligners within Snakefiles
"""


def hisat2_index_from_prefix(prefix):
    """
    Given a prefix, return a list of the corresponding hisat2 index files.
    """
    return ['{prefix}.{n}.ht2'.format(prefix=prefix, n=n) for n in range(1, 9)]


def prefix_from_hisat2_index(index_files):
    """
    Given a list of index files for hisat2, return the corresponding prefix.
    """
    if isinstance(index_files, str):
        return '.'.join(index_files.split('.')[:-2])
    else:
        prefixes = list(
            set(
                map(
                    lambda x: '.'.join(x.split('.')[:-2]), index_files)
            )
        )
        if len(prefixes) != 1:
            raise ValueError(
                "More than one prefix detected from '{0}'".format(prefixes)
            )
        return prefixes[0]


def bowtie2_index_from_prefix(prefix):
    """
    Given a prefix, return a list of the corresponding bowtie2 index files.
    """
    return (
        [
            '{prefix}.{n}.bt2'.format(prefix=prefix, n=n)
            for n in range(1, 5)
        ] + [
            '{prefix}.rev.{n}.bt2'.format(prefix=prefix, n=n)
            for n in range(1, 3)
        ]
    )


def prefix_from_bowtie2_index(index_files):
    """
    Given a list of index files for bowtie2, return the corresponding prefix.
    """
    if isinstance(index_files, str):
        return '.'.join(index_files.replace('.rev', '').split('.')[:-2])
    else:
        prefixes = list(
            set(
                map(
                    lambda x: '.'.join(x.replace('.rev', '').split('.')[:-2]),
                    index_files)
            )
        )
        if len(prefixes) != 1:
            raise ValueError(
                "More than one prefix detected from '{0}'".format(prefixes)
            )
        return prefixes[0]

def fastq_arg_from_input(fastqs):
    """
    Prepares the correct input FASTQ arguments for bowtie2 and HISAT2 based on
    whether or not the sample is paired-end.

    Parameters
    ----------
    fastqs : list-like
        List or snakemake.input object containing fastq filenames.
    """

    if isinstance(fastqs, str) or len(fastqs) == 1:
        fastqs = '-U {0} '.format(fastqs)
    else:
        assert len(fastqs) == 2
        fastqs = '-1 {0} -2 {1} '.format(*fastqs)
    return fastqs

