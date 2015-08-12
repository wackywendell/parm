from . import d2 as sim2
from . import d3 as sim3
from unittest import TestCase

from math import sqrt
import numpy as np

class VecTest(TestCase):
    def testAssignment(self):
        a = sim2.atom()
        v2 = a.x
        print('v2:', v2)
        x, y = v2
        self.assertAlmostEqual(x, 0.)
        self.assertAlmostEqual(y, 0.)
        a.x = (2.3, 3.6)
        x,y = a.x
        self.assertAlmostEqual(x, 2.3)
        self.assertAlmostEqual(y, 3.6)
    def testCross(self):
        x,y = sim2.cross((3.,-1.4), 2.)
        self.assertAlmostEqual(x, -2.8)
        self.assertAlmostEqual(y, -6.)
        c = sim2.cross((3.,-1.4), (-2., 0.))
        self.assertAlmostEqual(c, -2.8)
        c = sim2.cross((3.,-1.4), (0., -2.))
        self.assertAlmostEqual(c, -6.)
