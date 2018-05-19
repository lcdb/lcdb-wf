set -e
python -m doctest ../../ci/preprocessor.py
python ../../ci/preprocessor.py Snakefile > Snakefile.test && snakemake -s Snakefile.test "$@"
