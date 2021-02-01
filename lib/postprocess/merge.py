import os
from snakemake.shell import shell
from ..imports import resolve_name

def file_merge(origfns, newfn, *args):
    tmpfiles = ['{0}.{1}.sub.tmp'.format(newfn, i) for i in range(len(origfns))]
    try:
        for origfn, tmpfile, ppfunc in zip(origfns, tmpfiles, args):
            print(ppfunc)
            func = resolve_name(ppfunc)
            func(origfn, tmpfile)

        if os.path.exists(newfn):
            shell('rm {newfn}')

        if newfn.endswith('.gz'):
            fn = newfn.replace('.gz', '')
            for tmpfile in tmpfiles:
                shell("gunzip -c {tmpfile} >> {fn}")
            shell("gzip {fn}")
        else:
            for tmpfile in tmpfiles:
                shell("cat {tmpfile} >> {newfn}")

    except Exception as e:
        raise e

    finally:
        for i in tmpfiles:
            if os.path.exists(i):
                shell('rm {i}')

