import os
import pandas as pd
import pybedtools

README = """
The `peak_counts.tsv` file contains the number of called peaks in each of the
ChIP-seq peak-calling runs.
"""

outdir = os.path.dirname(snakemake.output[0])
if not os.path.exists(outdir):
    os.makedirs(outdir)

df = []
for i in snakemake.input:
    toks = i.split('/')
    peakcaller = toks[-3]
    label = toks[-2]
    df.append(dict(peakcaller=peakcaller, label=label, count=len(pybedtools.BedTool(i))))
(
    pd.DataFrame(df)[['peakcaller', 'label', 'count']]
    .to_csv(os.path.join(outdir, 'peak_counts.tsv'), sep='\t', index=False)
)

with open(snakemake.output[0], 'w') as fout:
    fout.write(README)
