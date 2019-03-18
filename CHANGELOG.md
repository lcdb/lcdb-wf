# v1.2
- updated requirements in `requirements.txt` and in wrappers
- changed all `pd.read_table()` to `pd.read_csv(sep="\t")` to prevent warnings
- changed all `yaml.load()` to `yaml.load(Loader=yaml.FullLoader)` to prevent warnings
- using DeprecationWarning rather than UserWarning in the deprecation handler
  so there's less spam in the logs
- from colocalization, removed the GAT "fractions" heatmap due to unresolved
  pandas index errors
- improved tests:
  - using data from pybedtools repo because modENCODE seems to be down
  - append rather than prepend base conda to PATH on circleci
