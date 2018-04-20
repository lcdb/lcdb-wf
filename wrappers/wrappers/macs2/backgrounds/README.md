# MACS2

Wraps several calls to the `macs2 bdgcmp` `macs2 bdgopt` `macs2 predictd` and `macs2 pileup` subprograms
for use in the middle of manual peak calling on bedgraph files. The benefit of this wrapper, as opposed
to callpeaks, is that (1) it supports somewhat more subtle modifications to macs2 internal parameters;
and (2) it caches bdg files according to which parameters are used to generate them, potentially
saving an extraordinary amount of time when you need to do parameter space scanning.