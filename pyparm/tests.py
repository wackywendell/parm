from . import d2 as sim2
from . import d3 as sim3
from unittest import TestCase
import textwrap

from math import sqrt
import numpy as np

class NPTestCase(TestCase):
    def assertClose(self, x, y, rtol=1e-05, atol=1e-08, msg=None):
        if np.allclose(x, y, rtol=rtol, atol=atol):
            return
        
        standardMsg = ('The following were not within %s rtol or %s atol:\n'
            '%s\n'
            '%s\n') % (rtol, atol, self._indent('A:  ', repr(x)), self._indent('B:  ', repr(y)))
        msg = self._formatMessage(msg, standardMsg)
        raise self.failureException(msg)
        
    def _indent(self, prefix, msg):
        lines = str(msg).splitlines()
        subsequent_indent = ' ' * len(prefix)
        lines =(
            [textwrap.wrap(l, width=80, 
                replace_whitespace=False,
                drop_whitespace=False,
                initial_indent=prefix,
                subsequent_indent=subsequent_indent,
                break_long_words=False,
                break_on_hyphens=False) for l in lines[:1]] +
            [textwrap.wrap(l, width=80, 
                replace_whitespace=False,
                drop_whitespace=False,
                initial_indent=subsequent_indent,
                subsequent_indent=subsequent_indent,
                break_long_words=False,
                break_on_hyphens=False) for l in lines[1:]])
        lines = [l for lset in lines for l in lset]
        return '\n'.join(lines)
        

class VecTest(NPTestCase):
    def test_assignment(self):
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
    def test_cross(self):
        x,y = sim2.cross((3.,-1.4), 2.)
        self.assertAlmostEqual(x, -2.8)
        self.assertAlmostEqual(y, -6.)
        c = sim2.cross((3.,-1.4), (-2., 0.))
        self.assertAlmostEqual(c, -2.8)
        c = sim2.cross((3.,-1.4), (0., -2.))
        self.assertAlmostEqual(c, -6.)

class RigidConstraintTest(NPTestCase):
    def setUp(self):
        self.masses8 = np.asarray([1.]*8)
        M = self.M8 = np.sum(self.masses8)
        dtheta = self.dtheta = np.pi/2.
        self.rot2d = np.asarray((
            [np.cos(dtheta), -np.sin(dtheta)],
            [np.sin(dtheta), np.cos(dtheta)]
        ))
        self.P8 = np.asarray(np.asarray([
            [0,0,0],
            [0,0,1],
            [0,1,0],
            [0,1,1],
            [1,0,0],
            [1,0,1],
            [1,1,0],
            [1,1,1.],
        ]))

        self.P8 -= np.sum(self.P8.T*self.masses8, axis=1).T/M

        x0, y0, z0 = self.P8.T
        x, y = self.rot2d.dot((x0, y0))

        self.Q8 = np.asarray((x, y, z0)).T
        
        self.infbox = sim3.InfiniteBox()
        self.atoms8 = sim3.atomvec(self.masses8)
        for a, loc in zip(self.atoms8, self.P8):
            a.x = loc
        self.rigid8 = sim3.RigidConstraint(self.infbox, self.atoms8)

        for a, loc in zip(self.atoms8, self.Q8):
            a.x = loc
    
    def test_cube_rot(self):
        m = self.rigid8.get_rotation()
        expected_m = np.array([
            [ 0., -1.,  0.],
            [ 1.,  0.,  0.],
            [ 0.,  0.,  1.]])
        self.assertClose(m, expected_m)


class RandomHertzianVerletTest(NPTestCase):
    phi = 0.3
    N_per_mol = 5
    N_mol = 12
    N = N_per_mol * N_mol
    dt = 0.01
    
    def setUp(self):
        # for consistency
        np.random.seed(131)
        self.radii = np.random.uniform(1.0, 2.0, size=(self.N,))
        self.masses = self.radii**3
        
        Vs = np.sum(self.radii**3)*4/3*np.pi
        V = Vs / self.phi
        self.L = float(V**(1./3.))
        box = self.box = sim3.OriginBox(self.L)
        
        self.atoms = sim3.atomvec(self.masses)
        self.hertz = sim3.Hertzian(self.box, self.atoms, 0.4)

        for a, radius, mass in zip(self.atoms, self.radii, self.masses):
            self.hertz.add(sim3.HertzianAtom(a, 100.0, radius*2.0, 2.0))
        
        self.reset_positions()
        
        collec = self.collec = sim3.collectionVerlet(self.box, self.atoms, self.dt, [self.hertz], [self.hertz.nlist()], [])
        
        print('EKUT:', collec.energy(), collec.kinetic(), collec.potentialenergy(), collec.temp())
    
    def reset_positions(self):
        np.random.seed(131)
        for a, radius, mass in zip(self.atoms, self.radii, self.masses):
            a.x = np.random.uniform(0., self.L, size=(3,))
            a.v = np.random.normal(size=(3,))
            a.f = np.random.normal(size=(3,))
        
        nl = self.hertz.nlist()
        nl.update_list(True)
    
    def reset(self):
        self.reset_positions()    
        self.collec.scaleVelocitiesT(1.0)
        for _ in range(1000):
            self.collec.timestep()
            self.collec.scaleVelocitiesT(1.0)
    
    def testEnergy(self):
        self.reset()
        collec = self.collec
        EKUTs = []
        lastE = collec.energy()
        for _ in range(1000):
            for _ in range(10):
                collec.timestep()
                self.assertClose(collec.energy(), lastE, rtol=1e-2)
                lastE = collec.energy()
            EKUTs.append((collec.energy(), collec.kinetic(), collec.potentialenergy(), collec.temp()))
            
        EKUTs = np.asarray(EKUTs)
        E,K,U,T = EKUTs.T
        
        mean_T = np.mean(T)
        self.assertClose(mean_T, 1.0, rtol=1e-1)
        
        som = np.std(E) / np.mean(E)
        self.assertClose(som, 0.0, atol=1e-2)


class RandomRigidConstraintTest(NPTestCase):
    phi_ish = 0.3
    N_per_mol = 4
    N_mol = 6
    N = N_per_mol * N_mol
    dt = 0.01
    
    def setUp(self):
        # for consistency
        np.random.seed(2460659162)
        self.radii = np.random.uniform(1.0, 2.0, size=(self.N,))
        self.masses = self.radii**3
        
        Vs = np.sum(self.radii**3)*4/3*np.pi
        V = Vs / self.phi_ish
        self.L = float(V**(1./3.))
        box = self.box = sim3.OriginBox(self.L)
        
        self.atoms = sim3.atomvec(self.masses)

        self.hertz = sim3.Hertzian(self.box, self.atoms, 0.4)

        self.subgroups = [sim3.subgroup(self.atoms) for _ in range(self.N_mol)]
        self.ixs = np.arange(self.N) // self.N_per_mol
        self.locs0 = np.random.normal(scale=10, size=(self.N_mol, 3))
        locs0 = np.array(self.locs0)
        for ix, a, radius, mass in zip(self.ixs, self.atoms, self.radii, self.masses):
            loc0 = locs0[ix]
            s = self.subgroups[ix]
            s.add(a)
            
            dx = np.random.normal(size=(3,))
            dx /= np.linalg.norm(dx)
            dx *= radius + 1.0
            a.x = loc0 + dx
            locs0[ix] = a.x
            a.v = np.random.normal(size=(3,))
            a.f = np.random.normal(size=(3,))
            self.hertz.add(sim3.HertzianAtom(a, 100.0, radius*2.0, 2.0))

        nl = self.hertz.nlist()
        for s in self.subgroups:
            for n,a in enumerate(s):
                alist = list(s)
                for m,a2 in enumerate(alist[:n]):
                    nl.ignore(a, a2)

        nl.update_list(True)
        
        self.rigids = [sim3.RigidConstraint(self.box, s) for s in self.subgroups]
        collec = self.collec = sim3.collectionVerlet(self.box, self.atoms, self.dt, [self.hertz], [self.hertz.nlist()], self.rigids)
        
        print('EKUT:', collec.energy(), collec.kinetic(), collec.potentialenergy(), collec.temp())
    
    def reset(self):
        np.random.seed(2460659162+1)
        locs0 = np.array(self.locs0)
        for ix, a, radius, mass in zip(self.ixs, self.atoms, self.radii, self.masses):
            loc0 = locs0[ix]
            s = self.subgroups[ix]
            
            dx = np.random.normal(size=(3,))
            dx /= np.linalg.norm(dx)
            dx *= radius + 1.0
            a.x = loc0 + dx
            locs0[ix] = a.x
            a.v = np.random.normal(size=(3,))
            a.f = np.random.normal(size=(3,))
        
        self.collec.scaleVelocitiesT(1.0)
        for _ in range(10000):
            self.collec.timestep()
            self.collec.scaleVelocitiesT(1.0)
    
    def testMomentOfInertia(self):
        self.reset()
        for s in self.subgroups:
            com = s.com()
            w = s.omega(com)
            L = s.angmomentum(com)
            I = s.moment(com)

            self.assertClose(L, I.dot(w))
    
    def testRigid(self):
        self.reset()
        collec = self.collec
        for _ in range(100):
            collec.timestep()
        for s, rigid in zip(self.subgroups, self.rigids):
            com = s.com()
            comv = s.comv()
            comf = s.comf()
            mom = s.moment(com)
            omega = s.omega(com)
            angm = s.angmomentum(com)
            torq = s.torque(com)
            K = s.kinetic(comv)
            KU = (collec.kinetic(), collec.potentialenergy())
            rigid.apply(self.box)

            com2 = s.com()
            comv2 = s.comv()
            comf2 = s.comf()
            angm2 = s.angmomentum(com2)
            mom2 = s.moment(com2)
            omega2 = s.omega(com2)
            torq2 = s.torque(com2)
            K2 = s.kinetic(comv2)
            KU2 = (collec.kinetic(), collec.potentialenergy())
            
            self.assertClose(com, com2)
            self.assertClose(comv, comv2)
            self.assertClose(comf, comf2)
            self.assertClose(K, K2)
            self.assertClose(KU, KU2)
            self.assertClose(comf, comf2)
            self.assertClose(angm, angm2)
            self.assertClose(mom, mom2)
            self.assertClose(omega, omega2)
            # Torque is actually not conserved.
            # self.assertClose(torq, torq2)
    
    def testEnergy(self):
        self.reset()
        collec = self.collec
        EKUTs = []
        lastE = collec.energy()
        for _ in range(1000):
            for _ in range(10):
                collec.timestep()
                self.assertClose(collec.energy(), lastE, rtol=1e-2)
                lastE = collec.energy()
            EKUTs.append((collec.energy(), collec.kinetic(), collec.potentialenergy(), collec.temp()))
            
        EKUTs = np.asarray(EKUTs)
        E,K,U,T = EKUTs.T
        
        mean_T = np.mean(T)
        self.assertClose(mean_T, 1.0, rtol=1e-1)
        
        som = np.std(E) / np.mean(E)
        self.assertClose(som, 0.0, atol=1e-2)
