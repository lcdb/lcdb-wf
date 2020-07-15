#!/usr/bin/env python

"""
Ensure that each named chunk in downstream Rmd files has a corresponding header
in the documentation.
"""


import sys
import os
import re

HERE = os.path.abspath(os.path.dirname(__file__))

# Keys are Rmd that have "```{r chunkname"; values are the corresponding .rst
# documentation.
d = {
    "../workflows/rnaseq/downstream/rnaseq.Rmd": "../docs/rnaseq-rmd.rst",
    "../workflows/rnaseq/downstream/functional-enrichment.Rmd": "../docs/functional-enrichment-rmarkdown-docs.rst",
    "../workflows/rnaseq/downstream/gene-patterns.Rmd": "../docs/gene-patterns-rmarkdown-docs.rst",
}

# For now, just assume "------" as H2 headings. A better approach would be to
# access the docutils parse tree
regex = re.compile(r"```\{r (?P<chunk>.*?)[\}, ]")


link_regex = re.compile(r"https://lcdb.github.io/lcdb-wf/rnaseq-rmd.html#(?P<name>)")

def get_chunk_names(rmd):
    lines = open(rmd).readlines()
    for i, line in enumerate(lines):
        m = regex.search(line)
        if not m:
            continue
        yield m["chunk"], lines[i+1].strip()


def get_headings(rst, underline="-"):
    last = None
    for line in open(rst):
        if set(line.strip()) == set("-"):
            yield last.replace('`', '')
        last = line.strip()

# H2 headings defined in this list are OK to not be present in the Rmd
ok_headings = ['Glossary']

errors = []

for rmd, rst in d.items():
    chunks, links = zip(*get_chunk_names(rmd))
    headings = set(get_headings(rst)).difference(ok_headings)
    chunks_without_headings = set(chunks).difference(headings)
    headings_without_chunks = headings.difference(chunks)
    if chunks_without_headings:
        errors.append(
            f"{rmd} has the following chunks undocumented in {rst}:"
            + "\n   - "
            + "\n   - ".join(chunks_without_headings)
            + "\n"
        )
    if headings_without_chunks:
        errors.append(
            f"{rst} has the following headings with no chunks in {rmd}:"
            + "\n   - "
            + "\n   - ".join(headings_without_chunks)
            + "\n"
        )

    for chunk, link in zip(chunks, links):
        if link_regex.search(link) is None:
            errors.append(
                f"In {rmd}, for chunk {chunk} the first line is:\n" +
                f"{link}\n" +
                f"Please change it to:\n" +
                f"# Docs: https://lcdb.github.io/lcdb-wf/rnaseq-rmd.html#{chunk}\n")


if errors:
    print("Identified chunks:")
    print(", ".join(sorted(chunks)))
    print("Identified headings:")
    print(", ".join(sorted(headings)))
    print("\n")
    print("\n".join(errors))
    sys.exit(1)
