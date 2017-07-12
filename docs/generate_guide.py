import os
from jinja2 import Template
HERE = os.path.abspath(os.path.dirname(__file__))
FILES = os.path.join(HERE, '..', 'guide-to-files.txt')

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

Guide to files
==============

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



"""#.format( '   '.join([i for i in open(FILES)]))


def setup(*args):
    t = Template(TEMPLATE)
    contents = t.render(files=files)
    with open('guide.rst', 'w') as fout:
        fout.write(contents)

if __name__ == "__main__":
    setup()
