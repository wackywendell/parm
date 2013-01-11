

import sys, math
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('files', nargs='+', type=str)
opts = parser.parse_args()

for fname in opts.files:
    ns = []

    with open(fname) as f:
        for l in f:
            if 'electricE' not in l:
                continue
            d = dict([s.split('=') for s in l.strip().split(' ')])
            ns.append(float(d['electricE']))
    
    if len(ns) == 0:
        print('%16s: (no data)' % fname)
        continue
    avg = sum(ns) / len(ns)
    std = math.sqrt((sum([n*n for n in ns]) / len(ns)) - avg*avg)

    print('%16s:%9.3f +-%9.3f' % (fname, avg, std))

