#!/usr/bin/env python2
from __future__ import print_function

import argparse, os, os.path, glob, shutil

parser = argparse.ArgumentParser()
parser.add_argument('files', nargs='+', type=str)
opts = parser.parse_args()

if len(opts.files) == 1 and os.path.isdir(opts.files[0]):
    d = opts.files[0]
    opts.files = glob.glob(d + '/*.xyz')


for fname in opts.files:
    bk = fname + '.bk'
    bkexists = os.path.exists(bk)
    readf = fname if not bkexists else bk
    print('reading', readf)
    with open(readf, 'r') as f:
        lines = f.readlines()

    print('finding')
    ls = [n for n,l in enumerate(lines) if 'time 0\n' in l]
    # go through and remove 'time 0' steps everywhere after the first one
    for n in reversed(ls[1:]):
        # remove 1015 lines starting at line n
        endn = n + 1015
        lines = lines[:n] + lines[endn:]
        print(n, endn)
        
    if len(ls) > 1:
        if os.path.exists(fname + '.bk'):
            print('backup exists')
        else:
            print('copying')
            shutil.copy2(fname, fname + '.bk')
        print('writing')
        with open(fname, 'w') as f:
            for l in lines:
                f.write(l)
    else:
        print('No issues found.')
        
