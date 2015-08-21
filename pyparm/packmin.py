import numpy as np

class Minimizer:
    def __init__(self, locs, diameters, masses=None, L=1.0, P=1e-4, dt=.1, CGerr=1e-12, Pfrac=1e-4,
                    need_contacts=False,
                    kappa=10.0, kmax=1000, secmax=40, 
                    seceps=1e-20, amax=2.0, dxmax=100, stepmax=1e-3,
                    itersteps=1000):
        """
        Minimizer to find a packing.
        
        Params
        ------
        itersteps : number of timesteps to take when using iter()
        CGerr : Maximum force magnitude allowed
        Pfrac : Allowed deviation from given pressure
        need_contacts : Require Nc >= Nc_exp to finish
        masses : mass of the particles; if None, will be diameters**3
        """
        self.itersteps = itersteps
        self.need_contacts = need_contacts
        self.CGerr = CGerr
        self.Pfrac = Pfrac
        self._L = L
        locs = np.array(locs)
        self._diameters = np.array(diameters)
        Ns, = self.diameters.shape
        Nl, ndim = np.shape(locs)
        if Nl != Ns:
            raise ValueError("Need shape N for diameters, Nx2 or Nx3 for locs; got {} and {}x{}".format(
                Ns, Nl, ndim))
        if ndim == 2:
            from . import d2 as sim
            self.sim = sim
        elif ndim == 3:
            from . import d3 as sim
            self.sim = sim
        else:
            raise ValueError("Number of dimensions must be 2 or 3; got {}".format(self.ndim))
        
        self.box = self.sim.OriginBox(L)
        
        self.masses = self.diameters**self.ndim if masses is None else masses
        self.atoms = self.sim.AtomVec([float(n) for n in self.masses])
        self.neighbors = self.sim.NeighborList(self.box, self.atoms, 0.4)
        self.hertz = self.sim.Hertzian(self.atoms, self.neighbors)

        for a, s, loc in zip(self.atoms, self.diameters, locs):
            a.x = self.sim.Vec(*loc)
            self.hertz.add(self.sim.HertzianAtom(a, 1.0, float(s), 2.0))
            
        collec = self.collec = self.sim.CollectionNLCG(self.box, self.atoms, dt, P, 
                [self.hertz], [self.neighbors], [], 
                kappa, kmax, secmax, seceps)
        collec.set_max_alpha(amax)
        collec.set_max_dx(dxmax)
        collec.set_max_step(stepmax)
        self.collec.set_forces(True, True)
        
        self.timesteps = 0
    
    @staticmethod
    def equal_mass(diameters, ndim):
        return [1] * len(diameters)
    
    @staticmethod
    def proportionate_mass(diameters, ndim):
        return np.array(diameters)**ndim
    
    @classmethod
    def randomized(cls, N=10, sizes=[1.0,1.4], ratios=None, ndim=3, phi0=0.01, 
                   mass_func=None, **kw):
        if mass_func is None:
            mass_func = cls.proportionate_mass
        if ratios == None:
            ratios = [1.] * len(sizes)
        if len(ratios) != len(sizes):
            raise ValueError("`ratios` list must be same length as `sizes` list")
        if ndim not in (2,3):
            raise ValueError("`ndim` must be 2 or 3")
        if not 0. < phi0 < 1.:
            raise ValueError("`phi0` must be 2 or 3")
        
        rmul = N / np.sum(ratios)

        ratios = np.array([r*rmul for r in ratios])
        rints, rfrac = np.array(np.floor(ratios), dtype=int), ratios % 1
        while np.sum(rints) < N:
            ix = np.argmax(rfrac)
            rints[ix] += 1
            rfrac[ix] = 0
        ratios = np.array(rints, dtype=int)
        assert sum(ratios) == N

        diameters = np.array([s for s,r in zip(sizes, ratios) for _ in range(r)])
        assert diameters.shape == (N,)
        
        masses = np.array(mass_func(diameters, ndim), dtype=float)
        assert masses.shape == (N,)
        
        Vs = np.sum(np.array(diameters)**ndim) * np.pi / (ndim*2)
        L = (Vs / phi0)**(1.0 / ndim)
        
        locs = np.random.rand(N, ndim) * L
        
        return cls(locs, diameters, masses, L=L, **kw)
        
    @property
    def diameters(self):
        return self._diameters
    
    @property
    def ndim(self):
        return self.sim.NDIM
    
    @property
    def L(self):
        return self.box.L()
    
    @L.setter
    def L(self, newL):
        self.box.resize_to_L(newL)
        self.collec.set_forces(True, True)
    
    def err(self):
        return (self.collec.pressure() / self.collec.P0 - 1, max([a.f.mag() for a in self.atoms]))
    
    def done(self):
        Perr, CGerr = self.err()
        errs = [abs(Perr) < self.Pfrac, CGerr < self.CGerr]
        if self.need_contacts:
            Nc, Nc_min, fl = self.pack_stats()
            errs.append(Nc > Nc_min and Nc_min > 4)
        
        return all(errs)
    
    def __iter__(self):
        while not self.done():
            yield next(self)
    
    def __next__(self):
        for _ in range(self.itersteps):
            self.collec.timestep()
        self.timesteps += self.itersteps
        return self.err()
    
    #-----------------------------------------------------------------------------------------------
    # Statistics (all intrinsic, i.e. divided by N)
    @property
    def N(self):
        return self.atoms.size()
    
    @property
    def Vspheres(self):
        return np.sum(self.diameters**self.sim.NDIM)*np.pi/(2*self.sim.NDIM*self.N)
    
    @property
    def V(self):
        return self.box.V() / self.N
    
    @property
    def phi(self):
        return self.Vspheres / self.V
    
    @property
    def H(self):
        return self.collec.hamiltonian() / self.N
    
    @property
    def U(self):
        return self.collec.potentialenergy() / self.N
    
    @property
    def pressure(self):
        return self.collec.pressure()
    
    @property
    def vdotv(self):
        return self.collec.vdotv() / self.N
    
    @property
    def overlap(self):
        "Average overlap per particle, in terms of distance"
        return self.hertz.energy(self.box) / self.N
    
    @property
    def locs(self):
        return np.array([a.x for a in self.atoms], dtype=float)
    
    def as_packing(self):
        from . import jammed
        L = self.L
        locs = np.remainder(self.locs + L/2., L) - L/2.
        
        return jammed.Packing(locs, self.diameters, L=L)
        
    def pack_stats(self):
        """Returns (number of backbone contacts, stable number, number of floaters)"""
        
        return self.as_packing().contacts()
    
    def status_str(self):
        Pdiff, CGerr = self.err()
        Nc, Nc_min, fl = self.pack_stats()
        Nc_min = max(Nc_min, 1)
        return 'dP: {:8.2g}, CGerr: {:7.2g}, phi: {:6.4f}. {:4d} Floaters, {:4d} / {:4d}'.format(
            Pdiff, CGerr, self.phi, fl, Nc, Nc_min)
