[![Build Status](https://travis-ci.org/lcdb/lcdb-wrapper-tests.svg?branch=master)](https://travis-ci.org/lcdb/lcdb-wrapper-tests)

# Testing Snakemake wrappers

This is a proof-of-concept for figuring out a good way of testing [Snakemake
wrappers](https://bitbucket.org/snakemake/snakemake-wrappers).

With the [`bioconda`](https://bioconda.github.io/) channel set up, install
the dependencies:

```
conda create -n testenv --file requirements.txt
```

And run tests:

```
source activate testenv
pytest test/tests.py -v
```

## Writing a test

Write new tests in `tests/test_<wrappername>.py`.

A test consists of a Snakefile (typically written as a triple-quoted string
inside each test function), an input data function (typically created via
`utils.symlink_in_tempdir`), and a function that performs arbitrary tests on
the output.

See `test/test_demo.py` for an annotated example.

## How the fixtures work

`py.test` *fixtures* are functions. They are used for things like setup and
teardown that needs to be done once, or for preparing data just once that can
be used across many tests. In this repo, the fixtures are primarily used for
downloading data from the testing data repo to a test-specific temporary
directory.

Fixtures are stored in a `conftest.py` file (see `tests/conftest.py`). Any
discovered pytest tests in any module (functions starting with `test_`) can use
any of these fixtures as an argument, and will get the return value of that
fixture at runtime.

pytest also has some built-in fixtures. Here we use the `tmpdir` fixture and
the `tmpdir_factory` fixture. The nice thing about these is that pytest
manages them on disk and only keeps the last handful of tempdirs. This makes it
straightforward to find the test data on disk and poke around when
troubleshooting.

Fixtures can use other fixtures, and what's really nice is that py.test handles
all of the dependencies across fixtures. For example, we can have a fixture
that downloads a FASTQ file and another fixture that runs FastQC on the FASTQ.
Then we can have a MultiQC test that uses the FastQC fixture as input. Py.test
will figure that out and create and clean up tempdirs as needed.

I'm still working out the best way to do this. As one extreme, we could
maintain a single big Snakefile for each set of mutually compatible wrappers,
and returns a big dict of all the files created. This does a high-level
functional test of all of the wrappers, while the more specific wrapper tests
could verify the contents.

The other extreme is to have a separate fixture for every stage of the
workflow.

### Specific example from the test suite

Let's take a concrete example. In the `tests/conftest.py` file, there's
a function called `sample1_se_fq`. Its job is to download a single FASTQ file
from the test data repo. It is marked with a `pytest.fixture` decorator, and is
set to run only once per session. Fixtures can use other fixtures, and
`sample1_se_fq()` uses the built-in `tmpdir_factory` fixture which is used for
once-per-session fixtures like we have here.

Once per test session (that is, a single invocation of the `pytest` command),
the `sample1_se_fq()` fixture will download the single-end FASTQ file to
`/tmp/pytest-of-$USER/pytest-$N/sample1_se_fq0/` (filling in the current user
and the current test number N). The fixture stores the path it downloaded and
can provide it to other fixtures or tests.

So that's the fixture. Then over in `tests/test_cutadapt.py`, we have tests for
the cutadapt wrapper. For example the `test_cutadapt_simple` test has two
fixtures: `sample1_se_fq` and `tmpdir`. So it will get the path to that
downloaded FASTQ returned by `sample1_se_fq` as its first argument, and will
get the built-in pytest fixture as its second argument. You can read the code
for details but basically this test copies the Snakefile and configured input
data (via the `symlink_in_tempdir` function) to the directory
`/tmp/pytest-of-$USER/pytest-$N/test_cutadapt_simple0/`.

**Importantly, this system of interdependent fixtures in temporary directories
means that for any test, you can go to
`/tmp/pytest-of-$USER/pytest-$N/$TEST_NAME` directory and you'll see the
Snakefile, input data, and any created output files.** This is useful for
writing the tests in the first place (e.g., what kinds of things to test for)
or for troubleshooting (e.g., go to the directory and run `snakemake` to see
what happens).

## Writing a wrapper test
First and foremost, read the source for examples; there are a fair number of
test cases now to use as templates.

When a new wrapper is added, add a new module to the `tests` dir named
`test_NAME.py` where NAME is the directory name under `wrappers`. That is, even
though there are wrappers for hisat2-build and hisat2-align, there's a single
`test_hisat2.py` module.

In that module, you'll probably want to start with these imports:

```python
from utils import run, dpath, symlink_in_tempdir
```

The test function should start with the name `test_`, and should at least use
the `tmpdir` fixture. It should use any other fixtures it might need. For
example, if we need a single-end FASTQ file, we can use the `sample1_se_fq`
fixture:

```python
def test_mywrapper(sample1_se_fq, tmpdir):
    pass
```

Now, a wrapper needs to be run within a Snakefile rule. And a rule generally
needs input data. So within a test function, write a Snakefile as a triple-quoted
string. You can use whatever input/output files you'd like, because the next
step will be setting those up.  When writing the rule, end it with the
following line:

```python
wrapper: 'file://wrapper'
```

This works because in a later step we'll be specifying which wrapper to use,
and it will be copied into a directory called `wrapper` alongside the Snakefile
in a tempdir.

Next we have to set up the input files. The input files will generally be
coming from a fixture. In the example described in the previous section, there
was a fixture that provided the filename of a FASTQ file downloaded to a temp
dir. We need to get that filename to the tempdir where the Snakefile and
wrapper are living. The `symlink_in_tempdir` does this. It takes a dictionary
as its only argument. This dict maps incoming fixture filenames to desired
symlinks. Use this to provide exactly the input files your Snakefile is expecting.

Technically, `symlink_in_tempdir` returns a function that takes a path as its
argument and symlinks keys over to values within that path. While this seems
a little convoluted, doing it this way means that we don't have to keep track
-- or even care -- what the fixture's provided filename is, avoiding the need
to keep looking back at the fixtures module to remember what the filenames are.
It keeps the input file setup logic tightly coupled to the Snakefile, since
they're both defined in the same function.

After mapping the input data, next we need to write a function. This function
will not get any arguments. It will be run in the same directory as the
Snakefile, right after the Snakefile is run, so you don't have to worry about
lots of `os.path.join` calls. A successful run of the Snakefile is one test,
but this is where you can test the actual contents of the output files.

Last, we tie everything together. Typically the last line of a test function looks like this:

```python
run(dpath('../wrappers/WRAPPERNAME`), snakefile, check, input_data_func, tmpdir)
```

`dpath` points to the wrapper as a path relative to the `tests` directory.

`snakefile` is the Snakefile as a string.

`check` is the function containing tests to run on the generated output.

`input_data_func` is the return value of `symlink_in_tempdir`.

`tmpdir` is the fixture coming in as one of this function's arguments.

The function is now discoverable by py.test. You can run pytest with the `-k`
option to select just this function to run. Once it passes, you should be good
to open a PR for testing on travis-ci.
