# -*- coding: utf-8 -*-

import sys, os.path
sys.path.append(os.path.expanduser('~/idp/src'))
import numpy as np
from math import sqrt

def rand_sphere(d0):
    """
    Get random points within a sphere of radius 1. Returns array of shape (d0, 3).
    """
    p1 = np.random.randn(d0, 3)
    m = np.sqrt(np.sum(p1**2, axis=1))
    
    rad = pow(np.random.rand(d0), 1.0/3.0)
    return (p1.T * (rad/m)).T

def rand_disk(d0):
    """
    Get random points within a disk of radius 1. Returns array of shape (d0, 2).
    """
    phi = np.random.rand(d0)*(2*np.pi)
    rad = pow(np.random.rand(d0), 1.0/2.0)
    return (np.array([cos(phi), sin(phi)]) * rad).T

def DM(rs, d):
    """Dynamical matrix for array rs, size ds. Assumes epsilon is the
    same for all.
    """
    N = len(d)
    rs = np.array(rs)
    M=np.zeros((2*N,2*N))
    for i in range(N):
        for j in range(i):
            rijvec=rs[i,:]-rs[j,:]
            rijvec=rijvec-np.around(rijvec)
            rijsq=np.sum(rijvec**2)
            dij=(d[i]+d[j])/2
            dijsq=dij**2
            if rijsq < dijsq:
                rij=np.sqrt(rijsq)
                rijouter = np.outer(rijvec,rijvec)
                Mij1=-rijouter/rijsq/dijsq
                Mij2=(1-rij/dij)*(rijouter/rijsq-np.eye(2))/rij/dij
                Mij=Mij1+Mij2
                
                M[2*i:2*i+2,2*j:2*j+2]=Mij
                M[2*j:2*j+2,2*i:2*i+2]=Mij
                M[2*i:2*i+2,2*i:2*i+2]-=Mij
                M[2*j:2*j+2,2*j:2*j+2]-=Mij
    return M

def DM_freqs(rs, ds):
    M = DM(rs,ds)
    (ew,ev) = np.linalg.eig(M)
    return np.sqrt(np.abs(ew)) / (2*np.pi)

def contacts(rs, ds, box):
    """Contact matrix for location Vecs rs, for particles of size ds in a box."""
    return np.array([
            [box.diff(r1,r2).mag() < (d1+d2)/2 for r1,d1 in zip(rs,ds)]
            for r2,d2 in zip(rs,ds)],
        dtype=bool)

def contact_change(cm1, cm2):
    """Returns number of contact (breaks, formations) given two contact maps."""
    return np.sum(cm1 & (~cm2)), np.sum(cm2 & (~cm1))

class Contacts:
    def __init__(self, atoms, sigmas, box):
        self.atoms = atoms
        self.sigmas = sigmas
        self.box = box
        rs = [a.x for a in self.atoms]
        self.cm = contacts(rs, self.sigmas, box=self.box)
    
    def update(self):
        old_cm = self.cm
        rs = [a.x for a in self.atoms]
        self.cm = contacts(rs, self.sigmas, box=self.box)
        return contact_change(old_cm, self.cm)
    
    def num_contacts(self):
        #~ print('num_contacts', np.sum(self.cm), '->', (np.sum(self.cm) - len(self.cm)) // 2)
        return (np.sum(self.cm) - len(self.cm)) // 2
    
def neighbors(xs,ys,sigmas, tol=1e-8):
    """
    For a set of particles at xs,ys with diameters sigmas, finds the 
    distance vector matrix (xdiff,ydiff) and the adjacency matrix.
    
    Assumes box size 1, returns (adjacency matrix, xdiff, ydiff)
    """
    xdiff = np.remainder(np.subtract.outer(xs, xs)+.5, 1)-.5
    ydiff = np.remainder(np.subtract.outer(ys, ys)+.5, 1)-.5
    sigmadists = np.add.outer(sigmas, sigmas)/2
    dists = np.sqrt((xdiff**2) + (ydiff**2))
    return dists - sigmadists < tol, xdiff, ydiff

def num_contacts(xs,ys,sigmas, tol=1e-8):
    matr, _, _ = neighbors(xs,ys,sigmas, tol=tol)
    N = len(sigmas)
    return (np.sum(matr) - N) // 2

def floater_indices(xs, ys, sigmas, tol=1e-8):
    areneighbors, _,_ = neighbors(xs,ys,sigmas,tol)
    notfloaters = np.sum(areneighbors, axis=0) >= 4 # 3 + itself
    #print("Floaters1", sum(notfloaters))
    oldNpack = -1
    Npack = np.sum(notfloaters)
    while Npack != oldNpack:
        areneighbors[~notfloaters] = 0
        areneighbors[:, ~notfloaters] = 0
        notfloaters = np.sum(areneighbors, axis=0) >= 4  # 3 + itself
        #print("Floaters2", sum(notfloaters))
        oldNpack, Npack = Npack, np.sum(notfloaters)
    return ~notfloaters

def PE(xs, ys, sigmas):
    N = len(sigmas)
    xdiffs = np.remainder(np.subtract.outer(xs, xs)+.5, 1)-.5
    ydiffs = np.remainder(np.subtract.outer(ys, ys)+.5, 1)-.5
    rdiffs = np.sqrt(xdiffs**2 + ydiffs**2)
    sigpairs = np.add.outer(sigmas,sigmas)/2.0
    
    # Only do i,j pairs for i<j
    rdiffs = rdiffs[np.triu_indices(N,1)]
    sigpairs = sigpairs[np.triu_indices(N,1)]
    
    # get max(sigma - dr, 0)
    drs = np.fmax(sigpairs - rdiffs, 0)
    
    return np.mean((drs/sigpairs)**2)/2.0

def P(xs, ys, sigmas):
    N = len(sigmas)
    xdiffs = np.remainder(np.subtract.outer(xs, xs)+.5, 1)-.5
    ydiffs = np.remainder(np.subtract.outer(ys, ys)+.5, 1)-.5
    
    rdiffs = np.sqrt(xdiffs**2 + ydiffs**2)
    sigpairs = np.add.outer(sigmas,sigmas)/2.0
    
    # Only do i,j pairs for i<j
    rdiffs = rdiffs[np.triu_indices(N,1)]
    sigpairs = sigpairs[np.triu_indices(N,1)]
    
    # get max(sigma - dr, 0)
    drs = np.fmax(sigpairs - rdiffs, 0)
    Fs = drs/sigpairs
    
    return np.mean(rdiffs*Fs)

def Fij_rij(xs, ys, sigmas):
    """Returns Fij, rij"""
    N = len(sigmas)
    didx = np.diag_indices(N)
    xdiffs = np.remainder(np.subtract.outer(xs, xs)+.5, 1)-.5
    ydiffs = np.remainder(np.subtract.outer(ys, ys)+.5, 1)-.5
    
    rij = np.array((xdiffs, ydiffs))
    rdiffs = np.sqrt(xdiffs**2 + ydiffs**2)
    sigpairs = np.add.outer(sigmas,sigmas)/2.0
    
    rhats = rij/rdiffs
    rhats[0][didx] = 0
    rhats[1][didx] = 0
    
    drs = np.fmax(1 - rdiffs/sigpairs, 0)
    drs[didx] = 0
    return rhats * (drs/sigpairs), rij

########################################################################
def neighbors3d(xs,ys,zs,sigmas, tol=1e-8):
    """
    Same as neighbors, but for 3D
    """
    xdiff = np.remainder(np.subtract.outer(xs, xs)+.5, 1)-.5
    ydiff = np.remainder(np.subtract.outer(ys, ys)+.5, 1)-.5
    zdiff = np.remainder(np.subtract.outer(zs, zs)+.5, 1)-.5
    sigmadists = np.add.outer(sigmas, sigmas)/2
    dists = np.sqrt((xdiff**2) + (ydiff**2) + (zdiff**2))
    return dists - sigmadists < tol, xdiff, ydiff, zdiff

def num_contacts3d(xs,ys,zs,sigmas, tol=1e-8):
    matr, _, _, _ = neighbors3d(xs,ys,zs,sigmas, tol=tol)
    N = len(sigmas)
    return (np.sum(matr) - N) // 2

def floater_indices3d(xs, ys, zs, sigmas, tol=1e-8):
    areneighbors, _,_,_ = neighbors3d(xs,ys,zs,sigmas,tol)
    notfloaters = np.sum(areneighbors, axis=0) >= 5
    
    oldNpack = -1
    Npack = np.sum(notfloaters)
    while Npack != oldNpack:
        areneighbors[~notfloaters] = 0
        areneighbors[:, ~notfloaters] = 0
        notfloaters = np.sum(areneighbors, axis=0) >= 5 # 4 + itself
        oldNpack, Npack = Npack, np.sum(notfloaters)
    
    return ~notfloaters

def nmer_neighbors3d(xs,ys,zs,sigmas, nmersize, tol=1e-8):
    """
    Same as neighbors, but for 3D
    """
    xdiff = np.remainder(np.subtract.outer(xs, xs)+.5, 1)-.5
    ydiff = np.remainder(np.subtract.outer(ys, ys)+.5, 1)-.5
    zdiff = np.remainder(np.subtract.outer(zs, zs)+.5, 1)-.5
    sigmadists = np.add.outer(sigmas, sigmas)/2
    dists = np.sqrt((xdiff**2) + (ydiff**2) + (zdiff**2))
    #~ print(dists, sigmadists)
    #~ exit()
    matr = dists - sigmadists < tol
    bigN = len(matr)
    N = bigN // int(nmersize)
    n = int(nmersize)
    smallmatr = np.zeros((N,N))
    for i in range(N):
        li,hi = i*n, (i+1)*n
        for j in range(N):
            lj,hj = j*n, (j+1)*n
            ijcontacts = np.sum(matr[li:hi, lj:hj])
            #~ print('ijcontacts:', i, j, matr[li:hi, lj:hj], ijcontacts)
            #~ print(np.shape(dists), np.shape(sigmadists))
            #~ print('dists - sigmadists:', dists[li:hi, lj:hj] - sigmadists[li:hi, lj:hj])
            if i == j: ijcontacts = 1
            smallmatr[i,j] = ijcontacts
    return smallmatr
            
def nmer_floater_indices3d(xs, ys, zs, sigmas, nmers, tol=1e-8):
    areneighbors = nmer_neighbors3d(xs,ys,zs,sigmas,nmers, tol=tol)
    contactlim = 6 if nmers >= 2 else 5 # 3 d.o.f. + (maybe) 2 rotational + 1 for one-sided
    notfloaters = np.sum(areneighbors, axis=0) >= contactlim + 1 # 1 more for itself
    
    oldNpack = -1
    Npack = np.sum(notfloaters)
    while Npack != oldNpack:
        areneighbors[~notfloaters] = 0
        areneighbors[:, ~notfloaters] = 0
        notfloaters = np.sum(areneighbors, axis=0) >= 7 # 6 + itself
        oldNpack, Npack = Npack, np.sum(notfloaters)
    
    return ~notfloaters

def nmer_contacts3d(xs,ys,zs,sigmas, nmers, tol=1e-8):
    matr = nmer_neighbors3d(xs,ys,zs,sigmas, nmers, tol=tol)
    N = len(matr)
    #~ print('contacts:',matr)
    return (np.sum(matr) - N) // 2

#-----------------------------------------------------------------------

def Vec2_diff(r1, r2, gamma=0.0, Lx = 1.0, Ly = 1.0):
    dx,dy = (r1 - r2).T
    im = np.round(dy/Ly)
    dy = dy - im*Ly
    dx = dx-np.round(dx/Lx-im*gamma)*Lx-im*gamma*Lx
    return np.array((dx,dy)).T

class Packing:
    def __init__(self, rs, sigmas, gamma = 0.0, L=1.0):
        self.rs = np.array(rs) / float(L)
        self.sigmas = np.array(sigmas) / float(L)
        self.gamma = gamma
        self.L = L
        
        self.N = len(self.sigmas)
        n, self.ndim = np.shape(self.rs)
        if n != self.N:
            raise ValueError("Need shape N for sigmas, Nx2 or Nx3 for rs; got {} and {}x{}".format(
                self.N, n, self.ndim))
        if self.ndim == 3:
            if self.gamma != 0:
                raise NotImplementedError("Can't do shearing for 3D")
        elif self.ndim != 2:
            raise ValueError("Number of dimensions must be 2 or 3; got {}".format(self.ndim))
    
    def neighbors(self, tol=1e-8):
        """
        For a set of particles at xs,ys with diameters sigmas, finds the 
        distance vector matrix (d x N x N) and the adjacency matrix.
        
        Assumes box size 1, returns (adjacency matrix, diffs)
        """
        diffs = np.remainder(np.array([np.subtract.outer(xs, xs) for xs in self.rs.T])+.5, 1) -.5
        
        if self.gamma != 0:
            xdiff, ydiff = diffs[:2]
            im = np.round(ydiff)
            xdiff -= im*self.gamma
            ydiff = ydiff - im
            xdiff -= np.round(xdiff)
            diffs[:2] = xdiff, ydiff
        
        sigmadists = np.add.outer(self.sigmas, self.sigmas)/2.
        dists = np.sqrt(np.sum(diffs**2, axis=0))
        
        return dists - sigmadists < tol, diffs*self.L
    
    def backbone(self, tol=1e-8):
        """Returns (backbone indices, neighbor matrix)"""
        areneighbors, _ = self.neighbors(tol)
        notfloaters = np.sum(areneighbors, axis=0) >= self.ndim + 2 # self.ndim + 1 for stability, +1 for itself
        
        oldNpack = -1
        Npack = np.sum(notfloaters)
        while Npack != oldNpack:
            areneighbors[~notfloaters] = 0
            areneighbors[:, ~notfloaters] = 0
            notfloaters = np.sum(areneighbors, axis=0) >= self.ndim + 2
            oldNpack, Npack = Npack, np.sum(notfloaters)
        
        return notfloaters, areneighbors
    
    def contacts(self, tol=1e-8):
        """Returns (number of backbone contacts, stable number, number of floaters)"""
        idx, nbor = self.backbone(tol=tol)
        return np.sum(np.triu(nbor, 1)), (np.sum(idx)-1) * self.ndim + 1, np.sum(~idx)
    
    def size_indices(self, tol=1e-8):
        """Returns [idx of sigma1, idx of sigma2, ...]"""
        sigs = np.array(np.round(self.sigmas/tol), dtype=int)
        sigset = set(sigs)
        return [sigs == s for s in sorted(sigset)]
    
    def dist_tree(self, other, tol=1e-8):
        if self.ndim == 2:
            from . import d2 as sim
        else:
            from . import d3 as sim
            
        assert self.ndim == other.ndim
        
        if self.gamma != 0 or other.gamma != 0:
            assert np.abs(self.gamma - other.gamma) <= tol
            gamma = (self.gamma + other.gamma) / 2.0
            box = sim.LeesEdwardsBox(sim.Vec(1,1), gamma)
        else:
            box = sim.OriginBox(1.0)
        
        sz1 = self.size_indices()
        assert(len(sz1) == 2)
        cutoff1 = int(np.sum(sz1[0]))
        
        sz2 = other.size_indices()
        assert(len(sz2) == 2)
        cutoff2 = int(np.sum(sz2[0]))
        
        vs1 = [sim.Vec(*xy) for idx in sz1 for xy in self.rs[idx]]
        vs2 = [sim.Vec(*xy) for idx in sz2 for xy in other.rs[idx]]
        
        tree = sim.JammingTreeBD(box, sim.vecvector(vs1), sim.vecvector(vs2), cutoff1, cutoff2)
        return tree
    
    def dist(self, other, tol=1e-8, maxt = 1000000):
        tree = self.dist_tree(other, tol=tol)
        tree.expand(maxt)
        return sqrt(tree.current_best().distance_squared)*self.L
        
    @staticmethod
    def _cage_pts(xyz, neighbor_xyzs, sigma, neighbor_sigmas, L, M, R):
        """Finds points within a distance R of point xyz that do not conflict with neigbors"""
        pts = rand_sphere(M)*R + xyz
        for nxyz, nsig in zip(neighbor_xyzs, neighbor_sigmas):
            dpts = np.remainder(pts - nxyz + L/2.0, L) - L/2.0
            dists_sq = np.sum(dpts**2, axis=1)
            goodix = dists_sq >= ((nsig + sigma)/2.0)**2
            pts = pts[goodix, :]
        return pts

    def cages(self, M=10000, R=None, Rfactor = 1.2, padding=0.1, Mfactor=0.1):
        """
        Find all cages in the current "packing". 
        
        The algorithm uses Monte Carlo: it finds M random points within a sphere of radius R from
        each particle, and sees if that particle could sit there without conflicting with other particles.
        Then (number of accepted points) / (number of test points) * (volume of sphere) is the
        volume of the cage.
        
        The algorithm is adaptive: if not enough test points are accepted (n < M * Mfactor), it tries
        more test points. If any test points are within `padding` of the edge, `R` is (temporarily)
        expanded.
        
        Parameters
        ----------
        M : Number of points in the sphere to test
        R : Size of sphere to test (will be expanded if necessary)
        Rfactor : How much to increase R by when the cage doesn't fit
        padding : How much larger the sphere should be than the cage (if it isn't, the sphere is 
                    expanded)
        Mfactor : Mfactor * M is the minimum number of points to find per cage. If they aren't 
                    found, more points are tested.
        
        Returns
        -------
        points : a list of (A x 3) lists, A indeterminate (but larger than M * Mfactor), with each 
                    list corresponding to the points within one cage.
        Vs : The approximate volumes of each cage.
        """
        if R is None: R = min(self.sigmas) * 0.2
        neighbordict = {}
        
        psets = []
        Vs = []
        for n, (xyz, s) in enumerate(zip(self.rs, self.sigmas)):
            curR = R
            curpow = -1
            nxyzs = nsigs = None
            
            def get_pts():
                pts = self._cage_pts(xyz, nxyzs, s, nsigs, 1.0, M, curR)
                max_dist = np.max(np.sqrt(np.sum((pts-xyz)**2, axis=1))) if len(pts) > 0 else 0
                return pts, max_dist
            
            pts, max_dist = [], curR
            while max_dist * (1. + padding) > curR:
                # print(n, curpow, max_dist * (1. + padding), curR)
                curpow += 1
                curR = R * pow(Rfactor, curpow)
                curM = M
                if curpow not in neighbordict:
                    pack = Packing(self.rs, self.sigmas + curR)
                    cur_neighbors, _ = pack.neighbors(tol=0)
                    cur_neighbors[np.diag_indices_from(cur_neighbors)] = False
                    neighbordict[curpow] = cur_neighbors
                cur_neighbors = neighbordict[curpow]
                nix = cur_neighbors[n]
                nxyzs = self.rs[nix, :]
                nsigs = self.sigmas[nix]
                pts, max_dist = get_pts()
                if max_dist * (1. + padding) > curR: continue
            
                while len(pts) < Mfactor * M:
                    # print(curM, len(pts), Mfactor * M, len(pts) < Mfactor * M)
                    pts2, maxdist2 = get_pts()
                    max_dist = max((max_dist, maxdist2))
                    pts = np.concatenate((pts, pts2))
                    curM += M
                    if max_dist * (1. + padding) > curR: break
                if max_dist > 0.5:
                    raise ValueError('Cage size filling entire space, cannot continue.')
            
            fracgood = len(pts) / curM
            r = fracgood**(1./3.) * curR
            V = fracgood * 4 * np.pi / 3 * (curR**3)
            psets.append(pts)
            Vs.append(V)
        return psets**self.L, np.array(Vs)*(self.L**self.ndim)
    
    def scene(pack, cmap=None, rot=0, camera_height = 0.7,  camera_dist = 1.5, angle=None, 
            lightstrength = 1.1, orthographic=False, pad=None, floatercolor=(.6,.6,.6)):
        """
        Render a scene.
        
        Requires `vapory` package, which requires the `povray` binary.
        
        Parameters
        ----------
        cmap : a colormap 
        
        Returns
        -------
        scene : vapory.Scene, which can be rendered using its `.render()` method.
        """
        import vapory
        import numpy as np
        
        try:
            import matplotlib as mpl
            import matplotlib.cm as mcm
            vmin, vmax = min(pack.sigmas), max(pack.sigmas)
            sm = mcm.ScalarMappable(norm=mpl.colors.Normalize(vmin, vmax), cmap=cmap)
            cols = [sm.to_rgba(s) for s in pack.sigmas]
        except ImportError:
            if not isinstance(cmap, list):
                raise ValueError("matplotlib could not be imported, and cmap not recognizeable as a list")
            cols = list(cmap)
        except TypeError:
            if not isinstance(cmap, list):
                raise ValueError("matplotlib could not convert cmap to a colormap, and cmap not recognizeable as a list")
            cols = list(cmap)

        if floatercolor is not None:
            ix, _ = pack.backbone()
            ns, = np.nonzero(~ix)
            for n in ns: cols[n] = floatercolor
        rs = np.remainder(pack.rs+.5, 1)-.5
        spheres = [
            vapory.Sphere(xyz, s/2., vapory.Texture( vapory.Pigment( 'color', col[:3] )))
            for xyz, s, col in zip(rs, pack.sigmas, cols)
        ]

        extent = (-.5, .5)
        corners = [np.array((x,y,z)) for x in extent for y in extent for z in extent]
        pairs = [(c1, c2) for c1 in corners for c2 in corners if np.allclose(np.sum((c1-c2)**2), 1) and sum(c1-c2) > 0]

        radius = 0.01
        col = vapory.Texture( vapory.Pigment( 'color', [.5,.5,.5]))
        cyls = [vapory.Cylinder(c1, c2, 0.01, col) for c1,c2 in pairs]
        caps = [vapory.Sphere(c, radius, col) for c in corners]
        
        
        light_locs = [
            [ 8.,  5., -3.],
            [-6.,  6., -5.],
            [-6., -7., -4.],
            [ 10., -5., 7.]
        ]
        rotlocs = [[x*np.cos(rot) - z*np.sin(rot),y, z*np.cos(rot) + x*np.sin(rot)] for x,y,z in light_locs]
        lights = [
            #vapory.LightSource( [2,3,5], 'color', [1,1,1] ),
            vapory.LightSource(loc, 'color', [lightstrength]*3 ) for loc in rotlocs
        ]
        cloc = [np.cos(rot)*camera_dist,camera_dist*camera_height,np.sin(rot)*camera_dist]
        mag = sqrt(sum([d**2 for d in cloc]))
        direction = [-v *2/ mag for v in cloc]
        
        if angle is None:
            if pad is None: pad = max(pack.sigmas)
            w = sqrt(2) + pad
            angle = float(np.arctan2(w, 2*camera_dist))*2*180/np.pi
        camera = vapory.Camera('location', cloc, 'look_at', [0,0,0], 'angle', angle)
        # vapory.Camera('orthographic', 'location', cloc, 'direction', direction, 'up', [0,2,0], 'right', [2,0,0])
         
        return vapory.Scene( camera, objects= lights + spheres + cyls + caps + [vapory.Background( "color", [1,1,1] )])
    
    def DM(self, masses=None):
        """Dynamical matrix for array rs, size ds. Assumes epsilon is the
        same for all.
        
        Parameters
        ----------
        masses : an array of length N of the masses of the particles.
        """
        N = len(self.sigmas)
        rs = self.rs
        d = self.ndim
        M=np.zeros((d*N,d*N))
        
        for i in range(N):
            sigi = self.sigmas[i]
            for j in range(i):
                rijvec=rs[i,:]-rs[j,:]
                rijvec=rijvec-np.around(rijvec)
                rijsq=np.sum(rijvec**2)
                dij=(sigi+self.sigmas[j])/2
                dijsq=dij**2
                if rijsq < dijsq:
                    rij=np.sqrt(rijsq)
                    rijouter = np.outer(rijvec,rijvec)
                    # U(r) = ½(1 - r/d)²	
                    # d²U/dxdy = (dr/dx)(dr/dy)/d² - (1 - r/d)(d²r/dxdy)/d
                    # dr/dx = x/r
                    # d²r/dxdy = -(x y) / r³
                    # d²U/dxdy = -(x y)/(r² d²) + (1 - r/d)((x y)/r²)/(d r)
                    # d²U/dx² = (dr/dx)²/d² - (1 - r/d)(d²r/dx²)/d
                    # d²r/dx² = -x² / r³ + 1/r
                    # d²U/dxᵢdxⱼ = -(xᵢ xⱼ)/(r² d²) + (1 - r/d)((xᵢ xⱼ)/r² - δᵢⱼ)/(d r)
                    
                    Mij1=-rijouter/rijsq/dijsq
                    Mij2=(1-rij/dij)*(rijouter/rijsq-np.eye(d))/rij/dij
                    Mij=Mij1+Mij2
                    
                    M[d*i:d*i+d,d*j:d*j+d]=Mij
                    M[d*j:d*j+d,d*i:d*i+d]=Mij
                    M[d*i:d*i+d,d*i:d*i+d]-=Mij
                    M[d*j:d*j+d,d*j:d*j+d]-=Mij
        
        np.divide(M, self.L**2, out=M)
        if masses is None:
            return M
        
        # TODO: is the mass part of this really part of this?
        marr = np.array(masses)
        assert np.shape(masses) == np.shape(self.sigmas)
        marr = np.array([masses]*d)
        marr = marr.T.flatten()
        # marr is now [m1,m1,m2,m2,...] (in 2D)
        mm = np.eye(d*N)
        np.multiply(mm, marr**-.5, out=mm)
        # mm is now M^-½, where M is the mass matrix
        
        mm.dot(M, out=M)
        M.dot(mm, out=M)
        return M
    
    def DM_freqs(self, masses=None):
        """Find the frequencies corresponding to the eigenvalues of the dynamical matrix.
        
        This is just a short wrapper around DM()."""
        ew, ev = np.linalg.eig(self.DM(masses=masses))
        # this used to be over 2pi; I don't know where the 2 went, but it seems to be gone now...
        return np.sqrt(np.abs(ew)) / (np.pi)

class Packing3d:
    def __init__(self, rs, sigmas, gamma = 0.0, L=1.0):
        self.rs = np.array(rs) / float(L)
        self.N, ndim = self.rs.shape
        assert(ndim == 3)
        self.sigmas = np.array(sigmas) / float(L)
    
    def neighbors(self, tol=1e-8):
        """
        The adjacency matrix.
        """
        xs,ys,zs = self.rs.T
        N = len(self.sigmas)
        arr = np.zeros((N,N), dtype=bool)
        for n,((x,y,z),s) in enumerate(zip(self.rs, self.sigmas)):
            xdiff = np.remainder(xs-x+.5, 1)-.5
            ydiff = np.remainder(ys-y+.5, 1)-.5
            zdiff = np.remainder(zs-z+.5, 1)-.5
            sigmadists = (self.sigmas + s)/2
            dists = np.sqrt((xdiff**2) + (ydiff**2) + (zdiff**2))
            arr[n] = dists - sigmadists < tol
        return arr
        
    def backbone(self, tol=1e-8):
        """Returns (backbone indices, neighbor matrix)"""
        areneighbors = self.neighbors(tol)
        notfloaters = np.sum(areneighbors, axis=0) >= 5 # 4 + itself
        
        oldNpack = -1
        Npack = np.sum(notfloaters)
        while Npack != oldNpack:
            areneighbors[~notfloaters] = 0
            areneighbors[:, ~notfloaters] = 0
            notfloaters = np.sum(areneighbors, axis=0) >= 5  # 4 + itself
            oldNpack, Npack = Npack, np.sum(notfloaters)
        
        return notfloaters, areneighbors
    
    def contacts(self, tol=1e-8):
        """Returns (number of backbone contacts, stable number)"""
        idx, nbor = self.backbone(tol=tol)
        return np.sum(np.triu(nbor)), np.sum(idx) * 3 - 2
    
    @staticmethod
    def _cage_pts(xyz, neighbor_xyzs, sigma, neighbor_sigmas, L, M, R):
        """Finds points within a distance R of point xyz that do not conflict with neigbors"""
        rpoints = rand_sphere(M) if self.ndim == 3 else rand_disk(M)
        pts = rpoints*R + xyz
        for nxyz, nsig in zip(neighbor_xyzs, neighbor_sigmas):
            dpts = np.remainder(pts - nxyz + L/2.0, L) - L/2.0
            dists_sq = np.sum(dpts**2, axis=1)
            goodix = dists_sq >= ((nsig + sigma)/2.0)**2
            pts = pts[goodix, :]
        return pts

    def cages(self, M=10000, R=None, Rfactor = 1.2, padding=0.1, Mfactor=0.1):
        """
        Find all cages in the current "packing". 
        
        The algorithm uses Monte Carlo: it finds M random points within a sphere of radius R from
        each particle, and sees if that particle could sit there without conflicting with other particles.
        Then (number of accepted points) / (number of test points) * (volume of sphere) is the
        volume of the cage.
        
        The algorithm is adaptive: if not enough test points are accepted (n < M * Mfactor), it tries
        more test points. If any test points are within `padding` of the edge, `R` is (temporarily)
        expanded.
        
        Parameters
        ----------
        M : Number of points in the sphere to test
        R : Size of sphere to test (will be expanded if necessary)
        Rfactor : How much to increase R by when the cage doesn't fit
        padding : How much larger the sphere should be than the cage (if it isn't, the sphere is 
                    expanded)
        Mfactor : Mfactor * M is the minimum number of points to find per cage. If they aren't 
                    found, more points are tested.
        
        Returns
        -------
        points : a list of (A x 3) lists, A indeterminate (but larger than M * Mfactor), with each 
                    list corresponding to the points within one cage.
        Vs : The approximate volumes of each cage.
        """
        if R is None: R = min(self.sigmas) * 0.2
        neighbordict = {}
        
        psets = []
        Vs = []
        for n, (xyz, s) in enumerate(zip(self.rs, self.sigmas)):
            curR = R
            curpow = -1
            nxyzs = nsigs = None
            
            def get_pts():
                pts = self._cage_pts(xyz, nxyzs, s, nsigs, 1.0, M, curR)
                max_dist = np.max(np.sqrt(np.sum((pts-xyz)**2, axis=1))) if len(pts) > 0 else 0
                return pts, max_dist
            
            pts, max_dist = [], curR
            while max_dist * (1. + padding) > curR:
                # print(n, curpow, max_dist * (1. + padding), curR)
                curpow += 1
                curR = R * pow(Rfactor, curpow)
                curM = M
                if curpow not in neighbordict:
                    pack = Packing3d(self.rs, self.sigmas + curR)
                    cur_neighbors = pack.neighbors3d(tol=0)
                    cur_neighbors[np.diag_indices_from(cur_neighbors)] = False
                    neighbordict[curpow] = cur_neighbors
                cur_neighbors = neighbordict[curpow]
                nix = cur_neighbors[n]
                nxyzs = self.rs[nix, :]
                nsigs = self.sigmas[nix]
                pts, max_dist = get_pts()
                if max_dist * (1. + padding) > curR: continue
            
                while len(pts) < Mfactor * M:
                    # print(curM, len(pts), Mfactor * M, len(pts) < Mfactor * M)
                    pts2, maxdist2 = get_pts()
                    max_dist = max((max_dist, maxdist2))
                    pts = np.concatenate((pts, pts2))
                    curM += M
                    if max_dist * (1. + padding) > curR: break
                if max_dist > 0.5:
                    raise ValueError('Cage size filling entire space, cannot continue.')
            
            fracgood = len(pts) / curM
            r = fracgood**(1./3.) * curR
            V = fracgood * 4 * np.pi / 3 * (curR**3)
            psets.append(pts)
            Vs.append(V)
        return psets, np.array(Vs)
