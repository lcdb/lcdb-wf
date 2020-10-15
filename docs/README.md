This documentation uses [sphinx](http://www.sphinx-doc.org) to buid the documentation.

The built documentation from the master branch can be found at
https://lcdb.github.io/lcdb-wf. If you want to build a local copy of the
documentation:

- create an environment from the `docs/docs-requirements.txt` file
- activate it
- run the Makefile in `docs`


That is:

```bash
# Create env
conda create -n lcdb-wf-docs \
  --file docs/docs-requirements.txt \
  --channel bioconda \
  --channel conda-forge \
  --channel lcdb

# activate it
source activate lcdb-wf-docs

# build the docs
cd docs
make html
```

The locally-built docs will be in `docs/_build/html/index.html`.
