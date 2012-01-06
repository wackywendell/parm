from __future__ import print_function

import cPickle as pickle
import os.path
import math

mydir = os.path.expanduser('~/idp/')
loadfile = mydir + 'blengths/stats.pkl'

rbonds, bbonds, rangles, bangles = pickle.load(open(loadfile))

bonds = [bbonds.values()] + [r.values() for r in rbonds.values()]
bonddicts = [bd for r in bonds for bd in r]
angles = [bangles.values()] + [r.values() for r in rangles.values()]
angledicts = [ad for r in angles for ad in r]

bondvariancesums = sum([b['std']**2 * b['N'] for b in bonddicts])
bondcount = sum([b['N'] for b in bonddicts])
bondstd = math.sqrt(bondvariancesums / bondcount)

print('bond std:', bondstd)
print('bond k: ', 1 / bondstd**2, 'T') 

angvariancesums = sum([b['std']**2 * b['N'] for b in angledicts])
angcount = sum([b['N'] for b in angledicts])
angstd = math.sqrt(angvariancesums / angcount)

print('angle std:', angstd)
print('angle k: ', 1 / angstd**2, 'T')
