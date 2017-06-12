#!/usr/bin/env python

"""

    Performs consistency checks for dependency versions.

    Checks requirements.txt against all the environment.yaml files from the
    wrappers. Prints out dependencies that are in wrappers but not defined in
    top-level requirements.txt, and prints out dependencies that are inconsistent
    among wrappers or between wrappers and top-level requirements.txt.

    Optionally writes a bash script with echo and sed commands that can be used to
    update everything to the latest version found collectively among requirements
    and wrappers.

    If a wrapper should be isolated from the top-level environment because its
    dependencies conflict (RSeQC I'm looking at you), then comment out those
    dependencies in requirements.txt rather than deleting them. The commands
    written the end will detect this, and will not try to add it to the list.

    However, if you bump a version in a commented-out line in requirements.txt,
    commands will be generated to bump those same versions in the wrappers.

    USE WITH CAUTION! This is trying to make your life easier by not having to
    dig through wrappers and update everything by hand, but definitely read
    through the proposed changes in the bash script to see what it wants to do
    and delete any offending lines.
"""

import yaml
import subprocess as sp
from collections import defaultdict
from distutils.version import LooseVersion


def toplevel_reqs(fn):
    """
    Returns a dictionary of {name: version} dependencies listed in `fn`.
    """
    reqs = {}
    for r in open('requirements.txt'):
        toks = r.lstrip('#').strip().split(' ')
        name = toks[0]
        if len(toks) == 1:
            version = 'undefined'
        else:
            version = toks[1]
        reqs[name] = {version: [fn]}
    return reqs


def wrapper_deps(d):
    """
    Returns a dictionary of {name: {version: [wrappers]}} found recursively in
    directory `d`.

    For example, this returned dictionary shows that the bowtie2 wrappers use
    a different version of bowtie2 than does fastq_screen while the RSeQC
    wrappers all agree on their pysam version::

        {
            'bowtie2': {
                '==2.2.6': ['fastq_screen'],
                '==2.3.0': ['bowtie2/align', 'bowtei2/build']
            },
            'cutadapt': {
                '==1.12': ['cutadapt']
            },
            'pysam': {
                '==0.10.0': ['geneBody_coverage', 'tin', 'bam_stat',
                'infer_experiment']
            },
        }
    """

    metas = (
        sp.check_output(
            ['find', d, '-name', 'environment.yaml'],
            universal_newlines=True)
        .splitlines(False)
    )
    rev = defaultdict(lambda: defaultdict(list))
    for fn in metas:
        m = yaml.load(open(fn))
        for dep in m['dependencies']:
            toks = dep.split(' ')
            name = toks[0]
            if len(toks) == 1:
                version = 'undefined'
            else:
                version = toks[1]
            rev[name][version].append(fn)
    for k, v in rev.items():
        rev[k] = dict(v)
    return rev

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(usage=__doc__)
    ap.add_argument('requirements', help='top-level requirements.txt')
    ap.add_argument('wrappers', help='wrappers directory')
    ap.add_argument('--cmds', help='Bash file to write out')

    args = ap.parse_args()
    rev = wrapper_deps(args.wrappers)
    reqs = toplevel_reqs(args.requirements)

    d = {}
    add_to_reqs = []
    keys = list(set(list(reqs.keys()) + list(rev.keys())))
    for k in keys:
        d[k] = defaultdict(list)
        if k not in reqs:
            add_to_reqs.append(' '.join([k] + list(rev[k].keys())))
        else:
            for k2, v2 in reqs[k].items():
                d[k][k2].extend(v2)
        if k not in rev:
            continue
        for k2, v2 in rev[k].items():
            d[k][k2].extend(v2)
    for k, v in d.items():
        d[k] = dict(v)

    cmds = ['#!/bin/bash']
    for i in sorted(add_to_reqs):
        cmds.append('echo "{}" >> requirements.txt'.format(i))
    cmds.append('sort requirements.txt > requirements.tmp && mv requirements.tmp requirements.txt')

    for k, v in d.items():
        if len(v) > 1:
            print('\n')
            recent = sorted(v.keys(), key=LooseVersion)[-1]
            for k2, v2 in v.items():
                print(k, k2)
                if k2 != recent:
                    for vi in v2:
                        cmds.append('sed -i "s/{k} {k2}/{k} {recent}/g" {vi}'.format(**locals()))
                print('\t' + '\n\t'.join(v2))

    if args.cmds:
        with open('cmds', 'w') as fout:
            fout.write('\n'.join(cmds))
    else:
        print('commands to run:')
        print()
        print('\n'.join(cmds))
