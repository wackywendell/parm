# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function, unicode_literals

from decimal import Decimal


def filefinder(dir, *names, regexp=None, matchall=True, types=Decimal,
               ext='tsv.gz', **args):
    """
    dir             : Directory to look in
    names           : keywords to look for / store into
    regexp          : Define the regular expression for groups
                      (default: #NAME##VALUE#-#NAME##VALUE#....ext)
    matchall        : Match the entire expression
    types           : What types to apply (default: all Decimal)
    ext             : File extension (ignores everything else) (default: 'tsv.gz')
    args            : Like names, but add some constants (such as T=1) for all fs
    """
    import re
    import fpath
    from collections import namedtuple
    
    names = [(n + '0' if n in 'f' else n) for n in names]
    FoundFile = namedtuple('FoundFile', names + ['f'])

    try:
        iter(types)
    except TypeError:
        types = [types] * len(names)
    
    dir = fpath.Dir(dir)
    children = [f for f in dir.children() if f[-1][(-len(ext)-1):] == '.' + ext]
    
    if regexp is None:
        regexp = '-'.join([n + r'(\F)' for n in names])
    regexp = regexp.replace(r'\F', r'(?:[0-9]*\.?[0-9]*)|(?:inf)')
    regexp = regexp.replace(r'\G', r'(?:[0-9]*\.?[0-9]*)|(?:inf)|(?:[0-9]*\.?[0-9]*e[-+0-9]+)')
    
    regex = re.compile(regexp)
    matches = [(list(regex.finditer(f[-1])), f) for f in children]
    badmatches = [f for ms, f in matches if len(ms) == 0]
    
    if badmatches and matchall:
        raise ValueError('%s could not match %s' % (regexp, badmatches[0][-1]))
    matches = [(ms[0], f) for ms, f in matches if len(ms) > 0]
    
    groups = sorted([list(zip(names, [t(m or 'nan') for t, m in
                     zip(types, mtch.groups())])) +
                     [('f', f)] for mtch, f in matches])
    
    pdicts = [dict(lst) for lst in groups]
    for p in pdicts:
        p.update(args)
    return [FoundFile(**d) for d in pdicts]


def filefinder2(dir, *names, delim='_', matchall=True, types=Decimal,
                ext='tsv.gz'):
    """
    Find files based on a standard naming scheme of name value [delim name value...] (literal .) ext
    
    dir             : Directory to look in
    names           : keywords to look for / store into
    delim           : Delimiter separating file names
    matchall        : Match every file with the given extension
    types           : What types to apply (default: all Decimal). May be a
                      single type or list of types
    ext             : File extension (ignores everything else) (default: 'tsv.gz')
    """
    import pathlib
    from collections import namedtuple
    FoundFile = namedtuple('FoundFile', list(names) + ['fname'])
    
    try:
        iter(types)
    except TypeError:
        types = [types] * len(names)
    
    dir = pathlib.Path(dir)
    children = [f for f in dir.iterdir() if f.name[(-len(ext)-1):] == '.' + ext]
    # print(len(children), [f[-1][(-len(ext)-1):] for f in dir.children()])
    
    def handlechild(c):
        cutf = c.name[:(-len(ext)-1)]
        splits = cutf.split(delim)
        if not len(splits) == len(names):
            raise ValueError("%r splits into %d instead of %d" %
                (cutf, len(splits), len(names)))
        ns, valstrs = zip(*[(v[:len(n)], v[len(n):]) for v, n in zip(splits, names)])
        m = 0
        
        nspace = dict()
        nspace['fname'] = c
        # nspace['cutf'] = cutf
        for n, v, n0, t in zip(ns, valstrs, names, types):
            m += 1
            if not n == n0:
                raise ValueError("%r does not match %r in %dth arg %r in %r" % (
                                 n, n0, m, n+v, cutf))
            try:
                value = t(v)
            except:
                raise ValueError("%r not converted by %r in %dth arg %r in %r" %
                                 (v, t, m, n+v, cutf))
            nspace[n0] = value
        return FoundFile(**nspace)
    
    def trycatch(c):
        try:
            return handlechild(c)
        except ValueError:
            return None
    
    if matchall:
        nspaces = [handlechild(c) for c in children]
    else:
        nspaces = [trycatch(c) for c in children]
        nspaces = [n for n in nspaces if n is not None]
        
    return sorted(nspaces)


def groupfiles(dlst, key):
    bigdict = dict()
    for d in dlst:
        curval = getattr(d, key)
        keylist = bigdict.get(curval, [])
        keylist.append(d)
        bigdict[curval] = keylist
    return bigdict


def groupby(dlst, groupkey, sortby=None):
    if sortby is None:
        return sorted(groupfiles(dlst, groupkey).items())
    else:
        return sorted(groupfiles(dlst, groupkey).items(), key=sortby)
