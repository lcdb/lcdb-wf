import os
import contextlib
import collections
from collections.abc import Iterable
from snakemake.shell import shell


@contextlib.contextmanager
def temp_env(env):
    """
    Context manager to temporarily set os.environ.
    """
    env = dict(env)
    orig = os.environ.copy()
    _env = {k: str(v) for k, v in env.items()}
    os.environ.update(_env)
    try:
        yield
    finally:
        os.environ.clear()
        os.environ.update(orig)


def flatten(iter, unlist=False):
    """
    Flatten an arbitrarily nested iterable whose innermost items are strings
    into a flat list of strings.

    Parameters
    ----------
    iter : iterable

    unlist : bool
        If True, convert single-item lists into a bare string
    """
    if isinstance(iter, dict):
        iter = iter.values()

    def gen():
        for item in iter:
            if isinstance(item, dict):
                item = item.values()
            if isinstance(item, Iterable) and not isinstance(item, str):
                yield from flatten(item)
            else:
                yield item
    results = list(gen())
    if unlist and len(results) == 1:
        return results[0]
    return results


def test_flatten():
    assert sorted(flatten({
        'a': {
            'b': {
                'c': ['a', 'b', 'c'],
            },
        },
        'x': ['e', 'f', 'g'],
        'y': {
            'z': 'd'
        },
    })) == ['a', 'b', 'c', 'd', 'e', 'f', 'g']

    assert flatten('a', True) == 'a'
    assert flatten(['a'], True) == 'a'
    assert flatten('a') == ['a']
    assert flatten(['a']) == ['a']


def updatecopy(orig, update_with, keys=None, override=False):
    """
    Update a copy of a dictionary, with a bit more control than the built-in
    dict.update.

    Parameters
    -----------

    orig : dict
        Dict to update

    update_with : dict
        Dict with new values

    keys : list or None
        If not None, then only consider these keys in `update_with`. Otherwise
        consider all.

    override : bool
        If True, then this is similar to `dict.update`, except only those keys
        in `keys` will be considered. If False (default), then if a key exists
        in both `orig` and `update_with`, no updating will occur so `orig` will
        retain its original value.
    """
    d = orig.copy()
    if keys is None:
        keys = update_with.keys()
    for k in keys:
        if k in update_with:
            if k in d and not override:
                continue
            d[k] = update_with[k]
    return d


def update_recursive(orig, update_with):
    """
    Recursively update one dict with another.

    From https://stackoverflow.com/a/3233356

    >>> orig = {'a': {'b': 1, 'c': 2, 'd': [7, 8, 9]}}
    >>> update_with = {'a': {'b': 5}}
    >>> expected = {'a': {'b': 5, 'c': 2, 'd': [7, 8, 9]}}
    >>> result = update_recursive(orig, update_with)
    >>> assert result == expected, result

    >>> update_with = {'a': {'d': 1}}
    >>> result = update_recursive(orig, update_with)
    >>> expected = {'a': {'b': 5, 'c': 2, 'd': 1}}
    >>> result = update_recursive(orig, update_with)
    >>> assert result == expected, result
    """
    for k, v in update_with.items():
        if isinstance(v, collections.Mapping):
            orig[k] = update_recursive(orig.get(k, {}), v)
        else:
            orig[k] = v
    return orig

def boolean_labels(names, idx, mapping={True: 'AND', False: 'NOT'},
                   strip='AND_'):
    """
    Creates labels for boolean lists.

    For example:

    >>> names = ['exp1', 'exp2', 'exp3']
    >>> idx = [True, True, False]
    >>> boolean_labels(names, idx)
    'exp1_AND_exp2_NOT_exp3'

    Parameters
    ----------

    names : list
        List of names to include in output

    idx : list
        List of booleans, same size as `names`

    mapping : dict
        Linking words to use for True and False

    strip : str
        Strip this text off the beginning of labels.

    given a list of names and a same-size boolean, return strings like

    a_NOT_b_AND_c

    or

    a_AND_b_AND_c_NOT_d_AND_e
    """
    s = []
    for i, (n, x) in enumerate(zip(names, idx)):
        s.append(mapping[x] + '_' + n)
    s = '_'.join(s)
    if s.startswith(strip):
        s = s.replace(strip, '', 1)
    return s


def make_relative_symlink(target, linkname):
    """
    Helper function to create a relative symlink.

    Changes to the dirname of the linkname and figures out the relative path to
    the target before creating the symlink.
    """
    linkdir = os.path.dirname(linkname)
    relative_target = os.path.relpath(target, start=linkdir)
    linkbase = os.path.basename(linkname)
    if not os.path.exists(linkdir):
        shell('mkdir -p {linkdir}')
    shell('cd {linkdir}; ln -sf {relative_target} {linkbase}')
