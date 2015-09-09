# -*- coding: utf-8 -*-

import sys
from datetime import datetime, timedelta


def norm(vec):
    from numpy.linalg import norm
    return norm(vec)


def sq(vec):
    from numpy.linalg import norm
    return norm(vec)**2


def printlist(lst, name, **kw):
    import numpy as np
    mean = float(np.mean(lst, 0))
    std = np.std(lst)
    stdstr = '%8.3f' % std if std < 10000 else '%8.5g' % std
    sovm = 100*std / mean if mean > 0 else float('nan')
    last = float(lst[-1])
    lstd = 100*(last - mean) / std if std > 0 else float('nan')
    print("{name:10s}={mean:10.5g}, "
          "σ={std:10.4g} ({sovm:9.3g}%) "
          "[{last:9.2f} ({lstd:7.3g}%)]".format(**locals()), **kw)


def to_dhms(tdelta):
    days = tdelta.days
    hr, sec = divmod(tdelta.seconds, 3600)
    mn, sec = divmod(sec, 60)
    return (days, hr, mn, sec)


def interval_str(days, hr, mn, sec):
    return (('%2dd ' % days) if days else '') + ('%2d:%02d:%02d' % (hr, mn, sec))


def roll(arr, axis):
    """'Roll' axis n to the front, so if it was A x B x C x D before, and
    axis was 3, its now C x A x B x D"""
    import numpy as np
    return np.rollaxis(arr, axis, 0)


def unroll(arr, axis):
    """Undo roll(arr, axis)"""
    import numpy as np
    arr = np.array(arr)
    axis = arr.ndim + axis + 1 if axis < 0 else axis+1
    return np.rollaxis(arr, 0, axis)
    
    
class Progress:
    def __init__(self, total):
        self.start = datetime.now()
        self.total = total
        
    def get_eta(self, i):
        """Returns (end time, time left) as datetime, timedelta objects."""
        now = datetime.now()
        tused = now - self.start
        
        # timedeltas cannot be multiplied by floats (?), so we get out the
        # seconds, use that, and make it back into a dt
        usedsecs = 24*60*60 * tused.days + tused.seconds + (1e-6 * tused.microseconds)
        tleft = (timedelta(0, (usedsecs) * ((self.total - i) / i)) if i > 0
                else timedelta(0, 0))
        
        return now + tleft, tleft
    
    def eta_strs(self, i, timefmt='(%a %b %d, %H:%M:%S)'):
        """Returns (time left, total time, endtime) as strs"""
        endtime, tleft = self.get_eta(i)
        interval = to_dhms(tleft)
        totinterval = to_dhms(endtime - self.start)
        
        return (
            interval_str(*interval),
            interval_str(*totinterval),
            endtime.strftime(timefmt)
        )
    
    def eta_str(self, i, timefmt='(%a %b %d, %H:%M:%S)'):
        ss = self.eta_strs(i, timefmt=timefmt)
        return 'Ends in ({} / {}) on {}'.format(*ss)


def errprint(*args, **kw):
    kwargs = {'file': sys.stderr}
    kwargs.update(kw)
    print(*args, **kwargs)
    kwargs['file'].flush()
