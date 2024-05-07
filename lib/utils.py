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
    assert (
        sorted(
            flatten(
                {
                    "a": {
                        "b": {
                            "c": ["a", "b", "c"],
                        },
                    },
                    "x": ["e", "f", "g"],
                    "y": {"z": "d"},
                }
            )
        )
        == ["a", "b", "c", "d", "e", "f", "g"]
    )

    assert flatten("a", True) == "a"
    assert flatten(["a"], True) == "a"
    assert flatten("a") == ["a"]
    assert flatten(["a"]) == ["a"]


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
        if isinstance(v, collections.abc.Mapping):
            orig[k] = update_recursive(orig.get(k, {}), v)
        else:
            orig[k] = v
    return orig


def boolean_labels(names, idx, mapping={True: "AND", False: "NOT"}, strip="AND_"):
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
        s.append(mapping[x] + "_" + n)
    s = "_".join(s)
    if s.startswith(strip):
        s = s.replace(strip, "", 1)
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
        shell("mkdir -p {linkdir}")
    shell("cd {linkdir}; ln -sf {relative_target} {linkbase}")


def autobump(*args, **kwargs):
    """
    Used to automatically bump resources depending on how many times the job
    was attempted. This will return a function that is appropriate to use for
    an entry in Snakemake's `resources:` directive::

        rule example:
            input: "a.txt"
            resources:
                mem_mb=autobump(gb=10),
                runtime=autobump(hours=2, increment_hours=10)

    Values can be specified in multiple ways.

    A single number will be provided as the resource, and will be used to
    increment each time. For example, this is the equivalent of 10 GB for the
    first attempt, and 20 GB for the second:

    >>> f = autobump(1024 * 10)
    >>> f(None, 1)
    10240

    Adding a second unnamed argument will use it as a value to increment by for
    each subsequent attempt. This will use 10 GB for the first attempt, and 110
    GB for the second attempt.

    >>> f = autobump(1024 * 10, 1024 * 100)
    >>> f(None, 1)
    10240

    >>> f(None, 2)
    112640

    Instead of bare numbers, keyword arguments can be used for more convenient
    specification of units. The above two examples can also take this form:

    >>> f = autobump(gb=10)
    >>> f(None, 1)
    10240

    >>> f = autobump(gb=10, increment_gb=100)
    >>> f(None, 2)
    112640


    Units can be minutes, hours, days, mb, gb, or tb. For example:

    >>> f = autobump(hours=2, increment_hours=5)
    >>> f(None, 2)
    420

    """
    multiplier = {
        "mb": 1,
        "minutes": 1,
        "gb": 1024,
        "hours": 60,
        "days": 1440,
        "tb": 1024 * 1024,
    }
    units = list(multiplier.keys())

    if args and kwargs:
        raise ValueError(
            "Mixture of unnamed and keyword arguments not supported with autobump()"
        )

    if len(kwargs) > 2:
        raise ValueError("Only 2 kwargs allowed for autobump()")

    elif len(args) == 1 and not kwargs:
        baseline_converted = args[0]
        increment_converted = baseline_converted

    elif len(args) == 2 and not kwargs:
        baseline_converted, increment_converted = args

    elif len(kwargs) <= 2:
        baseline_kwargs = [k for k in kwargs.keys() if k in units]
        if len(baseline_kwargs) != 1:
            raise ValueError(
                "Multiple baseline kwargs found. Do you need to change one to have an 'increment_' prefix?"
            )

        baseline_kwarg = baseline_kwargs[0]
        baseline_value = kwargs[baseline_kwarg]
        baseline_unit = baseline_kwarg

        increment_kwargs = [k for k in kwargs if k.startswith("increment_")]
        if increment_kwargs:
            assert len(increment_kwargs) == 1
            increment_kwarg = increment_kwargs[0]
            increment_value = kwargs[increment_kwarg]
            increment_unit = increment_kwarg.split("_")[-1]
        else:
            increment_value = baseline_value
            increment_unit = baseline_unit

        if baseline_unit not in multiplier:
            raise ValueError(
                f"Baseline unit {baseline_unit} not in valid units {units}"
            )
        if increment_unit not in multiplier:
            raise ValueError(
                f"Increment unit {increment_unit} not in valid units {units}"
            )

        baseline_converted = baseline_value * multiplier[baseline_unit]
        increment_converted = increment_value * multiplier[increment_unit]

    else:
        raise ValueError(f"Unhandled args and kwargs: {args}, {kwargs}")

    def f(wildcards, attempt):
        return  baseline_converted + (attempt - 1) * increment_converted

    return f


def gb(size_in_gb):
    return 1024 * size_in_gb


def hours(time_in_hours):
    return time_in_hours * 60
