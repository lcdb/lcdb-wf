set -e
python -m doctest ../../ci/preprocessor.py
../ci/preprocessor.py Snakefile > Snakefile.test && snakemake -s Snakefile.test "$@"
