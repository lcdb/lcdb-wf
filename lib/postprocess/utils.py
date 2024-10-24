import sys
import os
import re
import gzip
import zipfile
import shutil
import tempfile
import pandas as pd

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "../../lib"))
from common import openfile



def extract_from_zip(tmpfiles, outfile, path_in_zip):
    """
    Parameters
    ----------

    tmpfiles : list
        One-item list containing zip file

    outfile : str
        gzipped output file to create

    path_in_zip : str
        Path within zipfile to extract. You can identify the path using unzip
        -l x.zip from bash.
    """
    assert len(tmpfiles) == 1, f"expected single zip file, got {tmpfiles}"

    extraction_dir = tempfile.mkdtemp()

    with zipfile.ZipFile(tmpfiles[0], "r") as z:
        z.extract(path_in_zip, path=extraction_dir)

    full_path_to_extracted = os.path.join(extraction_dir, path_in_zip)

    with open(full_path_to_extracted, "rb") as fin:
        with gzip.open(outfile, "wb") as fout:
            shutil.copyfileobj(fin, fout)

    shutil.rmtree(extraction_dir)


def match_gtf_9th(tmpfiles, outfile, strmatch, optstrand = "None"):
    """
    Matches string to the 9th field of GTF and an optional strand that defaults to None;
    if the pattern is found and the provided strand match then the line is excluded

    Parameters
    ----------
    tmpfiles : str
        GTF files

    outfile : str
        gzipped output GTF file

    strmatch : list
        List of strings to match in the 9th field of the GTF. Must be list

    optstrand : str
        String to match to the strand. Default is None
    """
    regex_strmatch = re.compile(r'|'.join(strmatch))

    with gzip.open(outfile, 'wt') as fout:
        for tmpfn in tmpfiles:
            with openfile(tmpfn, 'rt') as tmp:
                for line in tmp:
                    if line.startswith("#"):
                        fout.write(line)
                    else:
                        toks = line.split('\t')
                        if not (regex_strmatch.search(toks[8]) != None and toks[6] == optstrand):
                            fout.write(line)

# match_gtf_9th(['/home/esnaultcm/Downloads/Rattus_norvegicus.Rnor_6.0.94.gtf.gz'], "test.gz", ['ENSRNOG00000046319'], '-')

def convert_gtf_chroms(tmpfiles, outfile, conv_table):
    """
    Convert chrom names in GTF file according to conversion table.

    Parameters
    ----------
    tmpfiles : str
        GTF files to look through

    outfile : str
        gzipped output GTF file

    conv_table : str
        Lookup table file for the chromosome name conversion. Uses pandas to
        read lookup table, so it can be file://, a path relative to the
        snakefile, or an http://, https://, or ftp:// URL.
    """
    lookup = pd.read_csv(
        conv_table, sep='\t', header=None, names=('a', 'b')
    ).set_index('a')['b'].to_dict()

    with gzip.open(outfile, 'wt') as fout:
        for tmpfn in tmpfiles:
            with openfile(tmpfn, 'rt') as tmp:
                for line in tmp:
                    if not line.startswith("#"):
                        toks = line.split('\t')
                        chrom = toks[0]
                        if chrom in lookup.keys():
                            toks[0]= lookup[chrom]
                            line = '\t'.join(toks)
                        else:
                            raise ValueError(
                                'Chromosome "{chrom}" not found in conversion table '
                                '"{conv_table}"'
                                .format(chrom=chrom, conv_table=conv_table)
                            )
                    fout.write(line)

def convert_fasta_chroms(tmpfiles, outfile, conv_table):
    """
    Convert chrom names in fasta file according to conversion table.

    Parameters
    ----------
    tmpfiles : str
        fasta files to look through

    outfile : str
        gzipped output fasta file

    conv_table : str
        Lookup table file for the chromosome name conversion. Uses pandas to
        read lookup table, so it can be file://, a path relative to the
        snakefile, or an http://, https://, or ftp:// URL.
    """

    lookup = pd.read_csv(
        conv_table, sep='\t', header=None, names=('a', 'b')
    ).set_index('a')['b'].to_dict()

    with gzip.open(outfile, 'wt') as fout:
        for tmpfn in tmpfiles:
            with openfile(tmpfn, 'rt') as tmp:
                for line in tmp:
                    if line.startswith(">"):
                        line = line.rstrip("\n")
                        toks = line.split(' ')
                        chrom = toks[0].lstrip(">")
                        chrom = chrom.rstrip("\n")
                        if chrom in lookup.keys():
                            toks[0]= ">" + lookup[chrom]
                            line = ' '.join(toks) + "\n"
                        else:
                            raise ValueError(
                                'Chromosome "{chrom}" not found in conversion table '
                                '"{conv_table}"'
                                .format(chrom=chrom, conv_table=conv_table)
                            )
                    fout.write(line)
