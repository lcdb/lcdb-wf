import gzip
import logging
import os
import re
import sys
import tempfile
import shutil
import zipfile

import gffutils
import pandas as pd
from snakemake.shell import shell

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "../../lib"))
from .. import utils as u

from . import *

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s", level=logging.INFO)


def default(origfn, newfn):
    shell("mv {origfn} {newfn}")


def ensure_single_unzipped(tmpfiles, outfile):
    """
    Sometimes it makes things easier in downstream code to assume there's
    a single uncompressed file to work with.
    """
    all_gzipped = all([u.is_gzipped(i) for i in tmpfiles])
    none_gzipped = all([not u.is_gzipped(i) for i in tmpfiles])

    if all_gzipped:
        shell("zcat {tmpfiles} > {outfile}")
        return outfile

    elif none_gzipped:
        shell("cat {tmpfiles} > {outfile}")
        return outfile

    else:
        raise ValueError("Mixture of compressed and uncompressed files")


def _patterns(include_patterns, exclude_patterns, verbose=False):
    """
    Return a function that will include/exclude strings based on the patterns
    provided.
    """

    if include_patterns and exclude_patterns:
        raise ValueError("include_patterns and exclude_patterns are mutually exclusive")
    patterns = []
    if include_patterns:
        for p in include_patterns:
            patterns.append(re.compile(p))

        def keep(s):     
            for p in patterns:
                if p.search(s):
                    if verbose:
                        logger.info(f"Keeping {s} because it matches {p}")
                    return True
            return False

    elif exclude_patterns:
        for p in exclude_patterns:
            patterns.append(re.compile(p))

        def keep(s):
            for p in patterns:
                if p.search(s):
                    if verbose:
                        logger.info(f"Excluding {s} because it matches {p}")
                    return False
            return True

    else:
        raise ValueError(
            "Expecting exactly one of include_patterns or exclude_patterns"
        )

    return keep


def filter_fasta_chroms(
    tmpfiles, outfile, include_patterns=None, exclude_patterns=None
):
    # samtools won't work with gzip (only bgzip) files, so the lowest common
    # denominator is to use uncompressed.
    working_file = ensure_single_unzipped(tmpfiles, outfile + ".tmp")
    if include_patterns and exclude_patterns:
        raise ValueError("include_patterns and exclude_patterns are mutually exclusive")

    logger.info(f"Finding chrom names and putting them in {working_file}.record_names")
    shell(
        'grep ">" {working_file} | cut -f1 -d " " | sed "s/>//g" > {working_file}.record_names'
    )

    keep = _patterns(include_patterns, exclude_patterns)
    with open(outfile + ".keep", "w") as fout, open(
        working_file + ".record_names", "r"
    ) as fin:
        for line in fin:
            line = line.replace(">", "").strip()
            chrom = line.split()[0]
            if keep(chrom):
                fout.write(chrom + "\n")
    shell("samtools faidx -r {outfile}.keep {working_file} | bgzip -c > {outfile}")
    # shell("rm {outfile}.tmp {outfile}.tmp.fai {outfile}.keep")
    shell("rm {tmpfiles}")


def filter_gtf_chroms(tmpfiles, outfile, include_patterns=None, exclude_patterns=None):
    working_file = ensure_single_unzipped(tmpfiles, outfile + ".tmp")
    keep = _patterns(include_patterns, exclude_patterns, verbose=False)
    with gzip.open(outfile, "wt") as fout:
        for feature in gffutils.DataIterator(working_file):
            if keep(feature.chrom):
                fout.write(str(feature) + "\n")
    shell("rm {tmpfiles}")


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


def match_gtf_9th(tmpfiles, outfile, strmatch, optstrand="None"):
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
    regex_strmatch = re.compile(r"|".join(strmatch))

    with gzip.open(outfile, "wt") as fout:
        for tmpfn in tmpfiles:
            with openfile(tmpfn, "rt") as tmp:
                for line in tmp:
                    if line.startswith("#"):
                        fout.write(line)
                    else:
                        toks = line.split("\t")
                        if not (
                            regex_strmatch.search(toks[8]) != None
                            and toks[6] == optstrand
                        ):
                            fout.write(line)



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
    lookup = (
        pd.read_csv(conv_table, sep="\t", header=None, names=("a", "b"))
        .set_index("a")["b"]
        .to_dict()
    )

    with gzip.open(outfile, "wt") as fout:
        for tmpfn in tmpfiles:
            with openfile(tmpfn, "rt") as tmp:
                for line in tmp:
                    if not line.startswith("#"):
                        toks = line.split("\t")
                        chrom = toks[0]
                        if chrom in lookup.keys():
                            toks[0] = lookup[chrom]
                            line = "\t".join(toks)
                        else:
                            raise ValueError(
                                'Chromosome "{chrom}" not found in conversion table '
                                '"{conv_table}"'.format(
                                    chrom=chrom, conv_table=conv_table
                                )
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

    lookup = (
        pd.read_csv(conv_table, sep="\t", header=None, names=("a", "b"))
        .set_index("a")["b"]
        .to_dict()
    )

    with gzip.open(outfile, "wt") as fout:
        for tmpfn in tmpfiles:
            with openfile(tmpfn, "rt") as tmp:
                for line in tmp:
                    if line.startswith(">"):
                        line = line.rstrip("\n")
                        toks = line.split(" ")
                        chrom = toks[0].lstrip(">")
                        chrom = chrom.rstrip("\n")
                        if chrom in lookup.keys():
                            toks[0] = ">" + lookup[chrom]
                            line = " ".join(toks) + "\n"
                        else:
                            raise ValueError(
                                'Chromosome "{chrom}" not found in conversion table '
                                '"{conv_table}"'.format(
                                    chrom=chrom, conv_table=conv_table
                                )
                            )
                    fout.write(line)
