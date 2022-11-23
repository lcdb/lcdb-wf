#!/bin/bash

set -e

# There are some settings in the Rmd files that only need to be used for
# testing, not production. These test-only lines are commented out with special
# strings. The preprocessor.py script edits the files to use the test settings.
# See the docstring of that file for details on how it works.
#
# Here, we run the preprocessor on all the Rmd files in downstream/, and store
# the newly-converted ones in dowstream-test/. Then we run rmarkdown::render on
# those new files.
mkdir -p downstream-test
for i in downstream/*.Rmd; do
    python ../../ci/preprocessor.py $i > downstream-test/$(basename $i)
done

# Make sure we move the config file there too
cp downstream/config.yaml downstream-test/config.yaml
cp downstream/text.yaml downstream-test/text.yaml

if [ -z "$ENV_R" ]; then
    echo Please set \$ENV_R to the conda env containing R and retry. Exiting.
    exit 1
fi

eval "$(conda shell.bash hook)"
conda activate $ENV_R
Rscript -e "rmarkdown::render('downstream-test/rnaseq.Rmd')"
conda deactivate
