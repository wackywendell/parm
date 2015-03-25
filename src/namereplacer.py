#!/bin/env python3
from pathlib import Path
import re
import sys

default_fs = (
        list(Path('src').glob('*.*pp')) + 
        list(Path('src/bin').glob('*.*pp')) + 
        [Path('src/sim.i')] +
        list(Path('pyparm').glob('*.py')) + 
        list(Path('pyparm/examples').glob('*.py')))

ignore_fs = [Path('pyparm') / 'd{}{}.py'.format(ndim, flt) for ndim in (2,3) for flt in ('', 'long')]

default_fs = [f for f in default_fs if f not in ignore_fs]

if len(sys.argv) > 1:
    # Passed command-line paths
    fs = []
    for p in sys.argv[1:]:
        p = Path(p)
        if not p.exists():
            print("Path {} not found".format(repr(p)))
            raise SystemExit
        fs.append(p)
else:
    fs = default_fs

with open('src/namereplacements.txt', 'r') as f:
    lines = [l.strip().replace('\t',' ').split(' ') for l in f]
    replacements = [(l[0], l[-1]) for l in lines if l[0] and l[-1]]

cnames = []
for f in fs:
    txt1 = txt = f.open('r').read()

    for old, new in replacements:
        txt = re.sub(r'\b\binteraction\b(?![.][hc]pp)\b'.format(old), new, txt)

    #newp = f.with_suffix(f.suffix + '.new')
    #print(newp)
    if txt1 != txt:
        print(f)
        with f.open('w') as newf:
            newf.write(txt)

    cnames.extend(re.findall(r'^\W*class ([a-z]\w*)', txt, re.M))
    cnames.extend(re.findall(r'^\W*struct ([a-z]\w*)', txt, re.M))

for cname in sorted(set(cnames)):
    print(cname, cname[0].capitalize() + cname[1:], sep='\t')
