# encoding: UTF-8
from sim import *

for name in list(locals().keys()):
    if '_swigregister' in name:
        del locals()[name]

class Stat:
    def __init__(self, name, func, shortname=None):
        self.name = name
        self.func = func
        self.shortname = shortname
        if not shortname:
            self.shortname = name
        self.values = []
        self.ts = []
    
    def update(self, t):
        self.ts.append(t)
        self.values.append(self.func())
    
    @staticmethod
    def plots(*stats):
        import numpy as np
        import matplotlib.pyplot as plt
        ymin = 0
        ymax = 0
        tmins = []
        tmaxs = []
        stdcol = '.75'
        labels = []
        valstart = .1
        for s in stats:
            startloc = int(len(s.values)*valstart)
            l = np.array(s.values[startloc:])
            ts = np.array(s.ts[startloc:])
            tmins.append(ts[0])
            tmaxs.append(ts[-1])
                
            mean = l.mean(0)
            std = l.std()
            sovm = 0
            if mean != 0:
                sovm = abs(std/mean)
            label = u'$\overline {{{}}} = {:.3f}$'.format(s.shortname,mean)
            stdlabel = '$\sigma={0:.3f}\;({1:.4f})$'.format(std, sovm)
            plt.axhline(mean, color=stdcol, linewidth=2)
            pstd=plt.axhline(mean-std, color=stdcol)
            plt.axhline(mean+std, color=stdcol)
            (p,)=plt.plot(ts, l, '.-')
            ymin = min(ymin, l.min())
            ymax = max(ymax, l.max())
            labels.append((p,label))
            labels.append((pstd,stdlabel))
        
        plt.axis([min(tmins), max(tmaxs), ymin, ymax ])
        plt.title(", ".join([s.name for s in stats]))
        plt.legend(*zip(*labels))
        plt.show()
    
    def plot(self):
        self.plots(self)
    
    def mean(self):
        import numpy as np
        return np.mean(self.values,0)
        
    def std(self):
        import numpy as np
        return np.std(self.values)
    
    def __str__(self):
        mean = self.mean()
        std = self.std()
        sovm = 0
        if mean != 0:
            sovm = abs(std*100.0/mean)
        return "{}={:.4g}+-{:.4g}%".format(self.name, mean, sovm)
