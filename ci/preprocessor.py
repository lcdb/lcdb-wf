#!/usr/bin/env python

"""
This script is used for preprocessing a file to prepare it for tests.

We often need to run tests with specific parameters that we may not want to use
in production. Rather than require users edit files to remove those
test-specific patterns, here we keep the test settings commented out and only
un-comment when running tests.

First, we look for any line that matches "# [test settings]" (case insensitive,
with optional surrounding spacing) and an optional signed integer. Any of these
would work:

    >>> assert matches('# [test settings]')
    >>> assert matches('#[test settings]')
    >>> assert matches('# [ test settings ]')
    >>> assert matches('# [ test settings -1]')
    >>> assert matches('# [ test settings +2]')
    >>> assert matches('# [ TEST SETTINGS +2]')
    >>> assert matches('# [ TeSt SeTTiNgS +2   ]')

If a lines does not match, output it as-is.

If a line matches, then uncomment it. Specifically, remove the first "#" in the
line; if it was followed by exactly one space, then remove that too.

If a line matches and a signed integer was provided, then consider it
a relative location, and then comment-out the referred-to line. Example:

    >>> preprocess('''
    ... use this for production
    ... # use this for tests  # [test settings -1]
    ... '''.splitlines(True))
    <BLANKLINE>
    # use this for production
    use this for tests  # [test settings -1]
    <BLANKLINE>

If the matched special string creates the first "#" in the line, then do
nothing to that line but still respect the relative locations. Useful for just
commenting out nearby lines for tests:

    >>> preprocess('''
    ... # [TEST SETTINGS +1]
    ... comment out for testing'''.splitlines(True))
    <BLANKLINE>
    # [TEST SETTINGS +1]
    # comment out for testing
"""

import re
regexp = re.compile(r'#\s?\[\s?test settings\s?(?P<rel>[-+]*\d)?\s*\]')


def matches(line):
    return regexp.search(line.lower()) is not None


def comment_line(line):
    """
    Adds a "#" just before the first non-whitespace character of a line.

    >>> assert comment_line('comment me') == '# comment me'
    >>> assert comment_line('    me too') == '    # me too'
    """
    x = []
    for i, character in enumerate(line):
        if character == ' ':
            x.append(character)
        else:
            break
    x.append('# ')
    x.extend(line[i:])
    return ''.join(x)


def uncomment_line(line):
    """
    Removes the first instance of "#" from a line; if it was followed by
    exactly one space then remove that too.

    >>> assert uncomment_line('# asdf') == 'asdf'
    >>> assert uncomment_line('#asdf') == 'asdf'
    >>> assert uncomment_line('# asdf # but this should be kept') == 'asdf # but this should be kept'
    >>> assert uncomment_line('#    asdf') == '    asdf'
    >>> assert uncomment_line('  #    asdf') == '      asdf'
    """
    first = line.find('#')

    # If the first comment is the one that flag the line, then do nothing.
    m = regexp.search(line.lower())
    if m:
        if m.start() == first:
            return line

    if line[first + 1] == ' ' and line[first + 2] != ' ':
        pattern = '# '
    else:
        pattern = '#'
    return line.replace(pattern, '', 1)


def preprocess(lines):

    if isinstance(lines, str):
        lines = [lines]

    # These lists will keep track of whether a line should be changed.  We need to
    # create them ahead of time so that we can use relative indexing from line N to
    # modify the state of lines N-1 or N+1
    uncomment = [False for i in range(len(lines))]
    comment = [False for i in range(len(lines))]

    for i, line in enumerate(lines):
        m = regexp.search(line.lower())
        if m:
            # There as at least a "[ test settings ]", so remove comment
            uncomment[i] = True

            # Figure out if there was also a relative location to uncomment,
            # and keep track of it in the `comment` list.
            rel = m.group('rel')
            if rel is not None:
                rel = int(rel)
                comment[i + rel] = True

    result = []
    for (c, u, line) in zip(comment, uncomment, lines):
        # E.g., in this situation, unclear what should happen:
        #
        #     # [test settings]
        #     # [test settings -1]
        #
        if c and u:
            raise ValueError("Line {0} is trying to be both commented and uncommented".format(line))
        if c:
            result.append(comment_line(line))
        elif u:
            result.append(uncomment_line(line))
        else:
            result.append(line)
    print(''.join(result))


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(usage=__doc__)
    ap.add_argument('infile', help='Input file to modify. Modified file printed to stdout.')
    args = ap.parse_args()
    lines = open(args.infile).readlines()
    preprocess(lines)
