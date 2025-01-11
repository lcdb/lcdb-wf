#!/usr/bin/env python

"""
This script is used for preprocessing a file to prepare it for tests.

We often need to run tests with specific parameters that we may not want to use
in production. Rather than require users edit files to remove those
test-specific patterns, here we keep the test settings commented out and only
un-comment when running tests.
"""


import re

regexp = re.compile(r"#\s?\[\s?(enable|disable) for test\s?\]")


def is_commented(line):
    return line.strip().startswith("#")


def comment_line(line):
    """
    Adds a "#" just before the first non-whitespace character of a line.

    >>> assert comment_line('comment me') == '# comment me'
    >>> assert comment_line('    me too') == '    # me too'
    """
    x = []
    for i, character in enumerate(line):
        if character == " ":
            x.append(character)
        else:
            break
    x.append("# ")
    x.extend(line[i:])
    return "".join(x)


def uncomment_line(line):
    """
    Removes the first instance of "#" from a line; if it was followed by
    exactly one space then remove that too . . . UNLESS the *only* comment is the
    special character that triggers this behavior, in which case we do nothing.

    >>> assert uncomment_line('# asdf') == 'asdf'
    >>> assert uncomment_line('#asdf') == 'asdf'
    >>> assert uncomment_line('# asdf # but this should be kept') == 'asdf # but this should be kept'
    >>> assert uncomment_line('#    asdf') == '    asdf'
    >>> assert uncomment_line('  #    asdf') == '      asdf'
    >>> assert uncomment_line('do nothing') == 'do nothing'
    >>> assert uncomment_line('do nothing # [disable for test]') == 'do nothing # [disable for test]')
    >>> assert uncomment_line('#uncomment # [disable for test]') == 'uncomment # [disable for test]')
    """
    first = line.find("#")

    # If the first comment is the one that flagged the line, then do nothing.
    m = regexp.search(line.lower())
    if m:
        if m.start() == first:
            return line

    if line[first + 1] == " " and line[first + 2] != " ":
        pattern = "# "
    else:
        pattern = "#"
    return line.replace(pattern, "", 1)


def preprocess(lines):
    result = []

    if isinstance(lines, str):
        lines = [lines]

    for line in lines:
        m = regexp.search(line.lower())
        if not m:
            result.append(line)
            continue

        action = m.group(1)
        if action == "enable" and is_commented(line):
            result.append(uncomment_line(line))
        elif action == "disable" and not is_commented(line):
            result.append(comment_line(line))
        else:
            raise ValueError(f"Inconsistent commenting and action:\n{line}")

    print("".join(result))


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(usage=__doc__)
    ap.add_argument(
        "infile", help="Input file to modify. Modified file printed to stdout."
    )
    args = ap.parse_args()
    lines = open(args.infile).readlines()
    preprocess(lines)
