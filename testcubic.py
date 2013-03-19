#!/usr/bin/python
from simw import solveCubic
import numpy as np

#~ a,b,c = np.random.normal(0,.1,(3,))
#~ poly = np.poly1d([1, a,b,c])
#~ x = solveCubic(a,b,c)
#~ 
#~ maxx = max(int(abs(x) + .5), 3)
#~ print('p(%.4f) = %.4g' % (x, poly(x)))
#~ rts = np.roots(poly)
#~ for r in rts:
    #~ if abs(r.imag) > 1e-8: continue
    #~ print('p(%.4f) = %.4g' % (r.real, poly(r.real)))
#~ 
#~ dx = min(abs(rts - x))
#~ if dx > 1e-8: print("BIG DIFFERENCE:", dx)

#~ from matplotlib import pyplot as plt
#~ xs = np.linspace(-maxx,maxx)
#~ plt.plot(xs, poly(xs), 'r-')
#~ plt.plot(x, poly(x), 'ko')
#~ plt.show()

mags = [1e-3, 1e-2, 1e-1, 1, 1e1, 1e2,1e3,1e4]
#~ mags = [1, 1e2]
for m1 in mags:
    for m2 in mags:
        #~ if m2 / m1 > 1e4 or m1 / m2 > 1e4: continue
        for m3 in mags:
          #~ if m3 / m2 > 1e4 or m2 / m3 > 1e4: continue
          for n in range(8):
            a = np.random.normal(0, m1)
            b = np.random.normal(0, m2)
            c = np.random.normal(0, m3)
            poly = np.poly1d([1, a,b,c])
            x = solveCubic(a,b,c)

            rts = np.roots(poly)
            
            dx = min(abs(rts - x))
            if dx > (1e-8 * np.max(np.abs([1,a,b,c]))):
                print("BIG DIFFERENCE:", dx)
                print("Magnitudes:", m1,m2,m3)
                print("Coeffs:    ", a,b,c)
                print("roots:     ", x, *rts)
                print("p(roots):  ", poly(x), *poly(rts))

a,b,c = 108.44300150958337, 118.30465621459336, 4.230956554439431
poly = np.poly1d([1, a,b,c])
x = solveCubic(a,b,c)
#~ print('p(x):', poly(x))
