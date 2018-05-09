python ../../ci/preprocessor.py Snakefile > Snakefile.test && snakemake -s Snakefile.test "$@"
