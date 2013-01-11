
import sys, xyzstats

fname = sys.argv[1]    
base, sep, ext = fname.rpartition('.')
statfname = base + sep + 'stats'
sk = xyzstats.statkeeper(fname, statfname)

import logging
#~ logging.getLogger().setLevel(logging.DEBUG)

Rgs = sk.Rg()
print('Found %d Rgs' % len(Rgs))

acs = sk.Rg_autocorr()

relax = sk.relax_Rg()
print('Relax:', relax)
if relax is not None: print('Times:', (sk.times[-1] - sk.times[0]) / relax)

from matplotlib import pyplot as plt

acs = acs[:(len(acs)/2)]
ts, vals = list(zip(*acs))

plt.plot(ts, vals)
#~ plt.ylim([-3,3])
plt.show()
