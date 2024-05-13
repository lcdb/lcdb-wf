#!/bin/bash

set -e

# There are some settings in the Rmd files that only need to be used for
# testing, not production. These test-only lines are commented out with special
# strings. The preprocessor.py script edits the files to use the test settings.
# See the docstring of that file for details on how it works.
#
# Here, we run the preprocessor on all the Rmd files in downstream/, and store
# the newly-converted ones in downstream-test/. Then we run rmarkdown::render on
# those new files.
TESTDIR=downstream-test
mkdir -p $TESTDIR
for i in downstream/*.Rmd; do
    python ../../ci/preprocessor.py $i > $TESTDIR/$(basename $i)
done

echo
echo "Standard workflow: rnaseq.Rmd -> functional-enrichment.Rmd -> gene-patterns.Rmd"

# Make sure we move the config file there too
cp downstream/config.yaml $TESTDIR/config.yaml
cp downstream/text.yaml $TESTDIR/text.yaml

# run rnaseq.Rmd
Rscript -e "rmarkdown::render('$TESTDIR/rnaseq.Rmd'); \
    sink('$TESTDIR/obj-names.txt'); cat(names(obj)); sink(); \
    sink('$TESTDIR/app-obj-names.txt'); cat(names(app_obj)); sink();"

OBS=$(cat $TESTDIR/obj-names.txt)
echo
echo "- check raw object slots: $OBS"
EXP="res_list dds_list"
[ "$OBS" != "$EXP" ] && echo "Object does not have expected slots" && exit 1

OBS=$(cat $TESTDIR/app-obj-names.txt)
echo "- check app object slots: $OBS"
EXP="res dds rld labels dds_mapping"
[ "$OBS" != "$EXP" ] && echo "Object does not have expected slots" && exit 1

# run functional-enrichment.Rmd
Rscript -e "rmarkdown::render('$TESTDIR/functional-enrichment.Rmd'); \
    sink('$TESTDIR/obj-names.txt'); cat(names(obj)); sink(); \
    sink('$TESTDIR/app-obj-names.txt'); cat(names(app_obj)); sink();"

OBS=$(cat $TESTDIR/obj-names.txt)
echo
echo "- check raw object slots: $OBS"
EXP="res_list dds_list rld_list enrich_list"
[ "$OBS" != "$EXP" ] && echo "Object does not have expected slots" && exit 1

OBS=$(cat $TESTDIR/app-obj-names.txt)
echo "- check app object slots: $OBS"
EXP="res dds rld labels dds_mapping enrich genetonic"
[ "$OBS" != "$EXP" ] && echo "Object does not have expected slots" && exit 1

# run gene-patterns.Rmd
Rscript -e "rmarkdown::render('$TESTDIR/gene-patterns.Rmd'); \
    sink('$TESTDIR/obj-names.txt'); cat(names(obj)); sink(); \
    sink('$TESTDIR/app-obj-names.txt'); cat(names(app_obj)); sink();"

OBS=$(cat $TESTDIR/obj-names.txt)
echo
echo "- check raw object slots: $OBS"
EXP="res_list dds_list rld_list degpatterns_list enrich_list"
[ "$OBS" != "$EXP" ] && echo "Object does not have expected slots" && exit 1

OBS=$(cat $TESTDIR/app-obj-names.txt)
echo "- check app object slots: $OBS"
EXP="res dds rld labels dds_mapping enrich genetonic degpatterns"
[ "$OBS" != "$EXP" ] && echo "Object does not have expected slots" && exit 1


echo
echo "Variation 1: rnaseq.Rmd -> gene-patterns.Rmd -> functional-enrichment.Rmd"

# now remove functional enrichment and gene patterns slots from combined-raw.Rds
# and rerun with those Rmds swapped
BEFORE=$(date -r $TESTDIR/final_clusters/consolidated_results.xlsx)

echo
echo "- Remove enrich_list & degpatterns_list slots from raw obj"
Rscript -e "obj <- readRDS('$TESTDIR/combined-raw.Rds'); \
    obj <- obj[ c('res_list', 'dds_list', 'rld_list') ]; \
    saveRDS(obj, '$TESTDIR/combined-raw.Rds', compress=FALSE)"

echo "- Remove enrich_list & degpatterns_list slots from app obj"
Rscript -e "obj <- readRDS('$TESTDIR/combined.Rds'); \
    obj <- obj[ c('res', 'dds', 'rld', 'labels', 'dds_mapping') ]; \
    saveRDS(obj, '$TESTDIR/combined.Rds', compress=FALSE)"

echo "- Empty functional enrichment and gene patterns cache"
rm -rf $TESTDIR/functional-enrichment_cache
rm -rf $TESTDIR/functional-enrichment_files
rm -rf $TESTDIR/gene-patterns_cache
rm -rf $TESTDIR/gene-patterns_files

Rscript -e "rmarkdown::render('$TESTDIR/gene-patterns.Rmd'); \
    sink('$TESTDIR/obj-names.txt'); cat(names(obj)); sink(); \
    sink('$TESTDIR/app-obj-names.txt'); cat(names(app_obj)); sink();"
OBS=$(cat $TESTDIR/obj-names.txt)
echo
echo "- check raw object slots: $OBS"
EXP="res_list dds_list rld_list degpatterns_list"
[ "$OBS" != "$EXP" ] && echo "Object does not have expected slots" && exit 1

OBS=$(cat $TESTDIR/app-obj-names.txt)
echo "- check app object slots: $OBS"
EXP="res dds rld labels dds_mapping degpatterns"
[ "$OBS" != "$EXP" ] && echo "Object does not have expected slots" && exit 1


Rscript -e "rmarkdown::render('$TESTDIR/functional-enrichment.Rmd'); \
    sink('$TESTDIR/obj-names.txt'); cat(names(obj)); sink(); \
    sink('$TESTDIR/app-obj-names.txt'); cat(names(app_obj)); sink();"
OBS=$(cat $TESTDIR/obj-names.txt)
echo
echo "- check raw object slots: $OBS"
EXP="res_list dds_list rld_list enrich_list degpatterns_list"
[ "$OBS" != "$EXP" ] && echo "Object does not have expected slots" && exit 1

OBS=$(cat $TESTDIR/app-obj-names.txt)
echo "- check app object slots: $OBS"
EXP="res dds rld labels dds_mapping enrich genetonic degpatterns"
[ "$OBS" != "$EXP" ] && echo "Object does not have expected slots" && exit 1

AFTER=$(date -r $TESTDIR/final_clusters/consolidated_results.xlsx)

echo
echo "- make sure results were updated"
[[ $AFTER < $BEFORE ]] && echo "- Results were not updated" && exit 1

# rerun after just touching raw RDS file
echo
echo "touch raw RDS file. Shouldn't update gene patterns"

BEFORE=$(date -r $TESTDIR/final_clusters/consolidated_results.xlsx)
touch $TESTDIR/combined-raw.Rds
Rscript -e "rmarkdown::render('$TESTDIR/gene-patterns.Rmd')"
AFTER=$(date -r $TESTDIR/final_clusters/consolidated_results.xlsx)
[[ $AFTER > $BEFORE ]] && echo "- Gene patterns was updated" && exit 1

echo "Done!"
