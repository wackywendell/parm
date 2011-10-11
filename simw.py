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
        stdcol = '.75'
        labels = []
        for s in stats:
            l = np.array(s.values)
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
            p=plt.plot(s.ts, l, '.-')
            ymin = min(ymin, l.min() * 1.1)
            ymax = max(ymax, 1.5*l.max())
            labels.append((p,label))
            labels.append((pstd,stdlabel))
        
        
        plt.axis([s.ts[-1]*.1, s.ts[-1], ymin, ymax ])
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
        return "{}={:.3f}+-{:.4f}%".format(self.name, mean, sovm)
