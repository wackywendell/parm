import numpy as np

class Minimizer:
    def __init__(self, locs, sigmas, masses=None, L=1.0, P=1e-4, dt=.1, CGerr=1e-12, Pfrac=1e-4,
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
        masses : mass of the particles; if None, will be sigmas**3
        """
        self.itersteps = itersteps
        self.need_contacts = need_contacts
        self.CGerr = CGerr
        self.Pfrac = Pfrac
        self._L = L
        locs = np.array(locs)
        self._sigmas = np.array(sigmas)
        Ns, = self.sigmas.shape
        Nl, ndim = np.shape(locs)
        if Nl != Ns:
            raise ValueError("Need shape N for sigmas, Nx2 or Nx3 for locs; got {} and {}x{}".format(
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
        
        self.masses = self.sigmas**self.ndim if masses is None else masses
        self.atoms = self.sim.atomvec([float(n) for n in self.masses])
        self.neighbors = self.sim.neighborlist(self.box, self.atoms, 0.4)
        self.hertz = self.sim.Hertzian(self.atoms, self.neighbors)

        for a, s, loc in zip(self.atoms, self.sigmas, locs):
            a.x = self.sim.Vec(*loc)
            self.hertz.add(self.sim.HertzianAtom(a, 1.0, float(s), 2.0))
            
        collec = self.collec = self.sim.collectionNLCG(self.box, self.atoms, dt, P, 
                [self.hertz], [self.neighbors], [], 
                kappa, kmax, secmax, seceps)
        collec.setamax(amax)
        collec.setdxmax(dxmax)
        collec.setstepmax(stepmax)
        self.collec.setForces(True, True)
        
        self.timesteps = 0
        
    @property
    def sigmas(self):
        return self._sigmas
    
    @property
    def ndim(self):
        return self.sim.NDIM
    
    @property
    def locs(self):
        return self._locs
    
    @locs.setter
    def locs(self, newlocs):
        newlocs = np.array(newlocs)
        assert locs.shape == newlocs.shape
        for a, loc in zip(self.atoms, newlocs):
            a.x = self.sim.Vec(*loc)
        self._locs = newlocs
        self.collec.setForces(True, True)
    
    @property
    def L(self):
        return self.box.L()
    
    @L.setter
    def L(self, newL):
        self.box.resizeL(newL)
        self.collec.setForces(True, True)
    
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
        return np.sum(self.sigmas**self.sim.NDIM)*np.pi/(2*self.sim.NDIM*self.N)
    
    @property
    def V(self):
        return self.box.V() / self.N
    
    @property
    def phi(self):
        return self.Vspheres / self.V
    
    @property
    def H(self):
        return self.collec.Hamiltonian() / self.N
    
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
        
        return jammed.Packing(locs, self.sigmas, L=L)
        
    def pack_stats(self):
        """Returns (number of backbone contacts, stable number, number of floaters)"""
        
        return self.as_packing().contacts()
    
    def status_str(self):
        Pdiff, CGerr = self.err()
        Nc, Nc_min, fl = self.pack_stats()
        Nc_min = max(Nc_min, 1)
        return 'dP: {:8.2g}, CGerr: {:7.2g}, phi: {:6.4f}. {:4d} Floaters, {:4d} / {:4d}'.format(
            Pdiff, CGerr, self.phi, fl, Nc, Nc_min)
