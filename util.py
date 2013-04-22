import sys
import numpy as np

def printlist(lst, name, **kw):
    mean = float(np.mean(lst,0))
    std = np.std(lst)
    stdstr = '%8.3f' % std if std < 10000 else '%8.5g' % std
    sovm = 100*std / mean if mean > 0 else float('nan')
    last = float(lst[-1])
    lstd = 100*(last - mean) / std if std > 0 else float('nan')
    print("{name:10s}={mean:10.5g}, σ={std:10.4g} ({sovm:9.3g}%) [{last:9.2f} ({lstd:7.3g}%)]".format(**locals()), **kw)

def to_dhms(tdelta):
    days = tdelta.days
    hr, sec = divmod(tdelta.seconds, 3600)
    mn, sec = divmod(sec, 60)
    return (days, hr, mn, sec)

def get_eta(cur, tot, starttime):
    from datetime import datetime, timedelta
    now = datetime.now()
    tused = now - starttime
    
    # timedeltas cannot be multiplied by floats (?), so we get out the 
    #seconds, use that, and make it back into a dt
    usedsecs = 24*60*60 * tused.days + tused.seconds
    tleft = timedelta(0, (usedsecs) * ((tot - cur) / cur))
    
    endtime = now + tleft
    return endtime, to_dhms(tleft)

def interval_str(days, hr, mn, sec):
    return (('%dd ' % days) if days else '') + ('%d:%02d:%02d'% (hr, mn, sec))

def errprint(*args, **kw):
    kwargs = {'file':sys.stderr}
    kwargs.update(kw)
    print(*args, **kwargs)
    kwargs['file'].flush()
