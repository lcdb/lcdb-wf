# This file demonstrates tests for the `demo` wrapper. It is heavily commented,
# and is included as part of the test suite to ensure that it's correct.

# The `run` function does most of the work. It creates a tempdir, copies over
# input data, Snakefile, and wrapper, runs the Snakefile, and runs
# a user-provided test function against the output.
from utils import run


# The `dpath` function figures out the path the wrapper even when in a tempdir
from utils import dpath

# `symlink_in_tempdir` is a decorator function that lets us easily map fixtures
# to input files expected by our Snakefile. The examples below will demonstrate
# how it works.
from utils import symlink_in_tempdir


# A note on fixtures
# ------------------
#
# py.test implicitly does a `from conftest import *`, so we will have the
# fixtures from that package available here.
#
# Currently we have the fixtures from raw_data_fixtures.py imported into
# conftest.py, which in turn makes them available in this file.
#
# py.test also includes a built-in `tmpdir` fixture which we use here to have
# a nicely-named tmpdir for running the test.
#
# See http://doc.pytest.org/en/latest/fixture.html for more info.


# Our first test. The test function names must start with `test_` in order for
# py.test to find them.
def test_demo(sample1_se_tiny_fq, tmpdir):

    # A note on these arguments
    # -------------------------
    #
    # Test function arguments are expected to be fixtures. The fixture
    # `sample1_se_tiny_fq` will be the path to the downloaded example data.  See
    # conftest.sample1_se_tiny_fq().
    #
    # The fixture `tmpdir` (which comes built-in with py.test) will be
    # a py.path.local object pointing to a tempdir created just for this test.
    # It will match the glob /tmp/pytest-*, and only the last 3 tempdirs are
    # retained.

    # Write the snakefile
    # -------------------
    # First we write the Snakefile to use in testing. Inputs need to come from
    # fixutres. Write whatever filename you'd like; we'll connect the fixture
    # to the written filename below.
    #
    # `snakefile` is typically a triple-quoted string; it will be automatically
    # run through textwrap.dedent later so you don't have to worry about
    # indentation.
    #
    # The wrapper will be copied to a subdirectory of the temp dir called,
    # appropriately enough, "wrapper". So your snakefile will generally end
    # with the line `wrapper: "file:wrapper"`.
    snakefile = '''
    rule demo:
        input: 'a.fastq.gz'
        output: 'b.fastq.gz'
        wrapper: "file:wrapper"
    '''

    # Map fixtures to input files
    # ---------------------------
    # Next we map the fixture sample1_se_tiny_fq (a temp file which has downloaded
    # data from the test data repo into a temp dir) to the input file that our
    # Snakefile expects.
    #
    # Keys are paths to downloaded example data (typically downloaded just once
    # per py.test session), which is provided by the fixture. The values of the
    # dict are paths relative to the Snakefile and must match what is expected
    # by the snakefile.
    #
    # Technically, `symlink_in_tempdir` returns a function that takes a path as
    # its argument and symlinks keys over to values within that path. While
    # this seems a little convoluted, doing it this way means that we don't
    # have to keep track -- or even care -- what the fixture's provided
    # filename is, avoiding the need to keep looking back at the fixtures
    # module to remember what the filenames are.  It keeps the input file setup
    # logic tightly coupled to the Snakefile, since they're both defined in the
    # same function.
    #
    # So: since the above snakefile expects a.fastq.gz as input, we need to
    # make that happen, like this:
    input_data_func=symlink_in_tempdir(
        {
            sample1_se_tiny_fq: 'a.fastq.gz'
        }
    )

    # Write a test function
    # ---------------------
    # This is our test function. It will be called after the Snakefile has been
    # run and it will be called in the same temp directory in which the
    # Snakefile is run, so paths should be relative to the Snakefile.
    #
    # This function should not accept any arguments.
    #
    # In this case, the demo wrapper simply copies input to output, so here we
    # assert the files are identical.
    def check():
        assert open('a.fastq.gz', 'rb').read() == open('b.fastq.gz', 'rb').read()

    # Call `run()`
    # ------------
    # Now that we have defined everything, the `run` function does all of the
    # work. Note we pass the `tmpdir` fixture here.
    #
    # (that's because py.test manages tmpdirs for tests, which are in this
    # current module, but run() lives in the utils module which won't get
    # nicely managed. But run() needs to know where to build the test case,
    # hence the need to pass it here)
    run(dpath('../wrappers/demo'), snakefile, check, input_data_func, tmpdir)



# This test function shows how to use downloaded paired-end data from
# a different fixture.
def test_demo_pe(sample1_pe_fq, tmpdir):

    # In contrast to the sample1_se_tiny_fq fixture used in the previous function,
    # here the paired-end fixture `sample1_pe_fq` is a tuple of path names (see
    # conftest.sample1_pe_fq())


    # The snakefile reflects what the wrapper expects for PE (see
    # wrappers/demo/README.md).
    snakefile = '''
    rule demo:
        input:
            R1='a1.fastq.gz',
            R2='a2.fastq.gz'
        output:
            R1='b1.fastq.gz',
            R2='b2.fastq.gz'
        wrapper: "file:wrapper"
    '''

    # Map fixture to input files. Again, since this is paired-end we need to
    # make sure both files are provided the right filename for testing.
    input_data_func=symlink_in_tempdir(
        {
            sample1_pe_fq[0]: 'a1.fastq.gz',
            sample1_pe_fq[1]: 'a2.fastq.gz',
        }
    )

    def check():
        assert open('a1.fastq.gz', 'rb').read() == open('b1.fastq.gz', 'rb').read()
        assert open('a2.fastq.gz', 'rb').read() == open('b2.fastq.gz', 'rb').read()

    run(dpath('../wrappers/demo'), snakefile, check, input_data_func, tmpdir)
