
import numpy as np
from random import random
from matplotlib import pylab as plt
from simw import autocorr as acorr

def norm(lst):
    l = np.array(lst)
    #~ l -= l.mean()
    l /= max(abs(l))
    return l

def aplot(lst,  title=None):
    plt.figure()
    plt.plot(lst)
    plt.plot(acorr(lst))
    if title:
        plt.title(title)
    plt.show()

def plotrand(n):
    m = np.cumsum(np.random.random(n) - .5)
    ac = acorr(m)
    m = norm(m)
    plt.figure()
    plt.plot(m)
    plt.plot(acorr(m))
    plt.title('random walk')
    plt.show()

plotrand(1000)
