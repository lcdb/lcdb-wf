import os
import pandas
import pybedtools
import gffutils
from gffutils import pybedtools_integration
db = gffutils.FeatureDB(snakemake.input.db)

outdir = os.path.dirname(snakemake.output[0])
if not os.path.exists(outdir):
    os.makedirs(outdir)

tsses = pybedtools_integration.tsses(db, as_bed6=True).saveas(os.path.join(outdir, 'tsses.bed'))

df = []
for pk in snakemake.input.peaks:
    name = os.path.basename(os.path.dirname(pk))
    p = pybedtools.BedTool(pk)
    peaks_at_tsses = p.intersect(tsses).saveas(os.path.join(outdir, name + '.bed'))
    npks = len(peaks_at_tsses)
    df.append(dict(name=name, npks=npks))

pandas.DataFrame(df)[['name', 'npks']].to_csv(snakemake.output[0], sep='\t', index=False)
