import os
from jinja2 import Template
HERE = os.path.abspath(os.path.dirname(__file__))
FILES = os.path.join(HERE, 'guide-to-files.txt')

class File(object):
    def __init__(self, fn):
        self._fn = fn.strip()
        self._desc = ""
        self._padding = 0

    @property
    def fn(self):
        f = self._fn
        if f.endswith('/'):
            cls = 'dir'
            padding = f.count('/') - 2
            f = os.path.basename(f.rstrip('/'))
            f = ('&nbsp;' * padding * 3) + f + '/'
        else:
            cls = 'file'
            d, f = os.path.split(f)
            padding = d.count('/')
            if self._fn.count('/') == 1:
                padding = 0
            f = ('&nbsp;' * padding * 3) + f

        if '***' in self._desc:
            cls += ' important'
            self._desc.replace('***', '')
        if not self.desc:
            cls += ' undoc'
        return '<span class="{0}">{1}</span>'.format(cls, f)

    @property
    def desc(self):
        return self._desc.replace('*', '')

    def __str__(self):
        if self.desc:
            return (
                '<a href="https://github.com/lcdb/lcdb-wf/blob/master/{0}" '
                'data-balloon="{1}" data-balloon-pos="right" '
                'data-balloon-length="xlarge" style="text-decoration:none;">{2}</a>'
                .format(self._fn, self.desc, self.fn)
            )
        return (
            '<a href="https://github.com/lcdb/lcdb-wf/blob/master/{0}" '
            'style="text-decoration:none;">{2}</a>'
            .format(self._fn, self.desc, self.fn)
        )


def gen():
    f = None
    for line in open(FILES):
        if line.startswith('/'):
            # it's a filename
            if f is not None:
                yield f
            f = File(line)
        else:
            f._desc += ' ' + line.strip()
    if f is not None:
        yield f

files = list(gen())


TEMPLATE = """\
.. _guide:

Guide to file hierarchy
=======================

The ``lcdb-wf`` workflow system is designed to have a standardized directory
structure and file hierarchy to allow us to be as consistent across many diverse
and disparate analyses and sources of data and reduce the overhead when it comes
to troubleshooting when something goes wrong. All the components of the repository
are laid out with this overarching design principle in mind.

Below we give a high-level overview and brief description of the files and folders used
by the workflows, and include an annotated directory tree highlighting the most important
parts of the repository.

Folder organization
~~~~~~~~~~~~~~~~~~~

The top level of the repo looks like this:

::

    [1]  ├── ci/
    [2]  ├── docs/
    [3]  ├── include/
    [4]  ├── lib/
    [5]  ├── README.md
    [6]  ├── requirements-non-r.txt
    [7]  ├── requirements-r.txt
    [8]  ├── workflows/
    [9]  └── wrappers/

1. ``ci`` contains infrastructure for continuous integration testing. You don't
   have to worry about this stuff unless you're actively developing `lcdb-wf`.

2. ``docs/`` contains the source for documentation. You're reading it.

3. ``include/`` has miscellaneous files and scripts that can be used by all
   workflows. Of particular note is the ``WRAPPER_SLURM`` script (see
   :ref:`cluster` for more) and the ``reference_configs`` directory (see
   :ref:`references` and :ref:`config` for more).

4. ``lib/`` contains Python modules used by the workflows.

5. ``README.md`` contains top-level info.

6. ``requirements-non-r.txt`` contains the package dependencies needed to run the
   workflows, and is used to set up a conda environment.

7. ``requirements-r.txt`` contains the package dependencies for R and various
   Bioconductor packages used in downstream analysis. See :ref:`conda-envs` for the
   rationale for splitting these.

8. ``workflows/`` contains one directory for each workflow. Each workflow directory contains
   its own ``Snakefile`` and configuration files. We go into more detail in the next section.

9. ``wrappers/`` contains Snakemake `wrappers
   <https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#wrappers>`_,
   which are scripts that can use their own independent environment. See
   :ref:`wrappers` for more.

Below, you can see a detailed overview of the files contained in these folders.


Annotated tree
~~~~~~~~~~~~~~

The following is an annotated directory tree of the ``lcdb-wf`` repository to
help orient you. Hover over files for a tooltip description; click a file to
view the most recent version on GitHub.

Files in bold are the most important.

.. raw:: html

    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/balloon-css/0.2.4/balloon.min.css">

    <style>
    .dir {
        font-family: monospace;
        font-size: 1em;
    }
    .file {
        font-family: monospace;
        font-size: 0.8em;
    }
    .important {
        font-weight: bold;
    }
    .undoc {
        color: #888;
        }

    </style>
    {% for x in files %}
    <p style="margin:0px;">{{ x }}</p>
    {% endfor %}


Now that you have seen which files and folders are the most important and have some idea
of where everything lives, let's look at how to run tests to make sure everything is set up 
correctly (see :ref:`running-the-tests`), or jump right in to learning about how to configure
the workflows for your particular experiment (see :ref:`config`).

"""#.format( '   '.join([i for i in open(FILES)]))


def setup(*args):
    t = Template(TEMPLATE)
    contents = t.render(files=files)
    with open('guide.rst', 'w') as fout:
        fout.write(contents)

if __name__ == "__main__":
    setup()
