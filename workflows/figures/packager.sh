#!/bin/bash
set -e
set -x
FILENAME=$(date +%Y%m%d_%H%M)_package.zip
OUTDIR=packages
mkdir -p $OUTDIR
rm -rf bundle
mkdir bundle

cat /dev/null > FILES
find figures >> FILES
find *html >> FILES
find *png >> FILES
rsync -rvL --files-from=FILES . bundle
zip -r "$FILENAME" bundle && mv "$FILENAME" "$OUTDIR" && rm -r bundle FILES
