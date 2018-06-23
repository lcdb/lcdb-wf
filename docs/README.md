This documentation uses [sphinx](http://www.sphinx-doc.org) to buid the documentation.

THe built documentation from the master branch can be found at
https://lcdb.github.io/lcdb-wf. If you want to build a local copy of the
documentation:

- create an environment
- activate it
- install sphinx into it
- run the Makefile in `docs`


That is:

```bash
# Create env
conda create -n lcdb-wf-docs \
  --file requirements.txt \
  --channel bioconda \
  --channel conda-forge \
  --channel lcdb

# activate it
source activate lcdb-wf-docs

# install sphinx
conda install sphinx

# build the docs
cd docs
make html
```

The locally-built docs will be in `docs/_build/html/index.html`.
