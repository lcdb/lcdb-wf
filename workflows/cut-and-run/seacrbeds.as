table mm10
"bigbed created from SEACR output"
(
string  chrom;		"Reference sequence chromosome or scaffold"
uint    chromStart;	"Start position of feature on chromosome"
uint    chromEnd;	"End position of feature on chromosome"
uint    totalSignal;		"Total signal contained within denoted coordinates"
uint    maxSignal;		"Maximum bedgraph signal attained at any base pair within denoted coordinates"
string  maxSignalLocation;		"Region representing the farthest upstream and farthest downstream bases within the denoted coordinates that are represented by the maximum bedgraph signal"
)
