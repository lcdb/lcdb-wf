table epic2InputPeak
"BED6+4 Peaks of signal enrichment based on pooled, normalized (interpreted) data."
(
    string chrom;        "Reference sequence chromosome or scaffold"
    uint   chromStart;   "Start position in chromosome"
    uint   chromEnd;     "End position in chromosome"
    float pValue;	 "Poisson-computed PValue based on the number of ChIP count vs. library-size normalized Input count in the region"
    uint   score;        "Indicates how dark the peak will be displayed in the browser (0-1000) "
    char[1]  strand;     "+ or - or . for unknown"
    float  name;  "The number of ChIP counts in the region (also including counts from windows with a count below the cutoff)"
    float  FDR;       "Benjamini-Hochberg correction of the p-values"
    float  log2FoldChange;       "Log2 of the region ChIP count vs. the library-size corrected region Input count"
)
