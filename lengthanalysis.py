
from Bio.PDB import *
from fpath import *
from collections import defaultdict

statfile = File('~/idp/blengths/stats.pkl').norm()
histfile = File('~/idp/blengths/hists.pkl').norm()

minbondlength = 2.5 # the largest mean is about 1.84, sigma=.02. 2.5 is way too much

bonddict = {
    'GLY' : [],
    'ALA' : [('CA','CB')],
    'SER' : [('CA','CB'), ('CB','OG')],
    'CYS' : [('CA','CB'), ('CB','SG')],
    'MET' : [('CA','CB'), ('CB','CG'), ('CG','SD'), ('SD','CE')],
    'LYS' : [('CA','CB'), ('CB','CG'), ('CG','CD'), ('CD','CE'), ('CE','NZ')],
    'VAL' : [('CA','CB'), ('CB','CG1'), ('CB','CG2')],
    'THR' : [('CA','CB'), ('CB','CG2'), ('CB','OG1')],
    'ILE' : [('CA','CB'), ('CB','CG1'), ('CB','CG2'), ('CG1','CD1')],
    'LEU' : [('CA','CB'), ('CB','CG'), ('CG','CD1'), ('CG','CD2')],
    'ASP' : [('CA','CB'), ('CB','CG'), ('CG','OD1'), ('CG','OD2')],
    'ASN' : [('CA','CB'), ('CB','CG'), ('CG','OD1'), ('CG','ND2')],
    'GLU' : [('CA','CB'), ('CB','CG'), ('CG','CD'), ('CD','OE1'), ('CD','OE2')],
    'GLN' : [('CA','CB'), ('CB','CG'), ('CG','CD'), ('CD','OE1'), ('CD','NE2')],
    'ARG' : [('CA','CB'), ('CB','CG'), ('CG','CD'), ('CD','NE'), 
                          ('NE','CZ'), ('CZ', 'NH1'), ('CZ', 'NH2')],
    'PRO' : [('CA','CB'), ('CB','CG'), ('CG','CD'), ('CD', 'N')],
    'HIS' : [('CA','CB'), ('CB','CG'), ('CG','ND1'), ('CG', 'CD2'),
                          ('ND1', 'CE1'), ('CD2', 'NE2'), ('CE1', 'NE2')],
    'PHE' : [('CA','CB'), ('CB','CG'), ('CG','CD1'), ('CG', 'CD2'),
             ('CD1', 'CE1'), ('CD2', 'CE2'), ('CE1', 'CZ'), ('CE2', 'CZ')],
    'TYR' : [('CA','CB'), ('CB','CG'), ('CG','CD1'), ('CG', 'CD2'),
             ('CD1', 'CE1'), ('CD2', 'CE2'), ('CE1', 'CZ'), 
             ('CE2', 'CZ'), ('CZ', 'OH')],
    'TRP' : [('CA','CB'), ('CB','CG'), ('CG','CD1'), ('CG', 'CD2'),
             ('CD1', 'NE1'), ('CD2', 'CE2'), ('CD2', 'CE3'), ('NE1', 'CE2'), 
             ('CE2', 'CZ2'), ('CE3', 'CZ3'), ('CZ2', 'CH2'), ('CZ3', 'CH2')]
    }

otherbonds = [('N', 'CA'), ('CA', 'C'), ('C', 'O')]
lastbond = ('C', 'OXT')
pairbond  = ('C', 'N')

def get_bonds(resname, atom):
    allbonds = bonddict[resname] + list(otherbonds) + [lastbond]
    for (a1,a2) in allbonds:
        if a1 == atom:
            yield a2
        elif a2 == atom:
            yield a1

# atoms per residue
atom_names = defaultdict(set)
for res,bonds in list(bonddict.items()):
    for a1, a2 in bonds:
        atom_names[res].add(a1)
        atom_names[res].add(a2)
    atom_names[res].add('N')
    atom_names[res].add('CA')
    atom_names[res].add('C')
    atom_names[res].add('O')

# all possible bond angles per residue
angles = dict()
for res, atomset in list(atom_names.items()):
    curset = set()
    for atom in atomset:
        bonds = list(get_bonds(res, atom))
        # for every pair of atoms bonded to this atom, add it as an angle
        for b1 in bonds:
            for b2 in bonds:
                if b2 > b1:
                    curset.add((b1,atom,b2))
    angles[res] = curset

pairangles = (('C','N','CA'), # last this this
              ('O','C','N'), # this this next
              ('CA','C','N') # this this next
              )
lastangle = ('CA','C','OXT')

def getangle(res, a1, a2, a3):
    try:
        vs = [res[a].get_vector() for a in (a1,a2,a3)]
    except KeyError:
        return None
    return calc_angle(*vs)

broken = []

def getlengths(res, lastres=None):
    err = False # whether any atoms were missing
    if(res.resname not in bonddict):
        raise KeyError(res.resname + ' NOT VALID RESNAME')
    resdict = dict()
    backdict = dict()
    
    def distFindAppend(a1, a2, d, thiserr):
        try:
            dist = res[a1] - res[a2]
        except KeyError:
            if thiserr:
                broken.append((res, a1, a2))
                err = True
        else:
            if dist > minbondlength:
                err = True
                return
            d[(a1, a2)] = dist
    
    # internal bonds
    for a1, a2 in bonddict[res.resname]:
        distFindAppend(a1, a2, resdict, True)
    
    #backbone bonds
    for a1, a2 in otherbonds:
        distFindAppend(a1, a2, backdict, True)
    
    #possible chain-ending bond
    a1, a2 = lastbond
    distFindAppend(a1, a2, backdict, False)
    
    # neighbor-neighbor peptide bond
    if lastres:
        a1,a2 = pairbond
        try:
            dist = lastres[a1] - res[a2]
        except KeyError:
            broken.append((res, a1, a2))
            err = True
        else:
            if dist < minbondlength: backdict[pairbond] = dist

    return resdict, backdict, err

def getangles(res, lastres=None, nextres=None):
    err = False # whether any atoms were missing
    if(res.resname not in angles):
        raise KeyError(res.resname + ' NOT VALID RESNAME')
    resdict = dict()
    backdict = dict()
    
    # internal and backbone bonds
    #~ print('---')
    for bondatoms in angles[res.resname]:
        ang = getangle(res, *bondatoms)
        angdict = resdict
        if all(a in ('C','O','CA','N','OXT') for a in bondatoms):
            #~ if(ang is not None): print("IN:", *bondatoms)
            angdict = backdict
        if ang is not None:
            angdict[bondatoms] = ang
        else:
            err = True
    
    #possible last angle
    ang = getangle(res, *lastangle)
    if ang is not None:
        backdict[lastangle] = ang
    
    #peptide-bond angles
    def pairang(atoms, resA, resB):
        if resA is None or resB is None: return
        a1, a2, a3 = atoms
        try: a1, a2, a3 = resA[a1], res[a2], resB[a3]
        except KeyError: err=True
        else:
            #~ print(a1,a2,a3)
            if a2 - a1 < minbondlength and a3 - a2 < minbondlength:
                vs = [a.get_vector() for a in (a1,a2,a3)]
                backdict[atoms] = calc_angle(*vs)
                #~ print(a1,a2,a3,backdict[pairangles[0]])
            else:
                print(resA, resB, a2-a1, a3-a2, *atoms)
                err=True
    
    pairang(pairangles[0], lastres, res)
    pairang(pairangles[1], res, nextres)
    pairang(pairangles[2], res, nextres)
    
    #~ if lastres is not None:
        #~ a1, a2, a3 = pairangles[0]
        #~ try: a1, a2, a3 = lastres[a1], res[a2], res[a3]
        #~ except KeyError: err=True
        #~ else:
            #~ if a2 - a1 < 10:
                #~ vs = [a.get_vector() for a in (a1,a2,a3)]
                #~ backdict[pairangles[0]] = calc_angle(*vs)
            #~ else: err=True
    #~ if nextres is not None:
        #~ a1, a2, a3 = pairangles[1]
        #~ try: a1, a2, a3 = res[a1], res[a2], nextres[a3]
        #~ except KeyError: err=True
        #~ else:
            #~ if a3 - a2 < 10:
                #~ vs = [a.get_vector() for a in (a1,a2,a3)]
                #~ backdict[pairangles[1]] = calc_angle(*vs)
            #~ else: err = True
    
    return resdict, backdict, err


pdbp = PDBParser(QUIET=True)
backbone_bonds = defaultdict(list)
res_bonds = defaultdict(lambda: defaultdict(list))
res_angles = defaultdict(lambda: defaultdict(list))
back_angles = defaultdict(list)
notfound = defaultdict(int)

def res_dist(res1, res2):
    hetflag, resseq1, icode=res1.get_id()
    hetflag, resseq2, icode=res2.get_id()
    return resseq1 - resseq2

def parsepdb(fname, permissive):
    global bdict, resdict
    struc = pdbp.get_structure('temp',fname)
    chain = list(struc.get_chains())[0]
    lastres = None
    reslist = [None] + chain.child_list + [None]
    lastfail, curfail = False, False
    for rs in zip(reslist, reslist[1:-1], reslist[2:]):
        lastres, res, nextres = rs
        if res.resname not in bonddict:
            notfound[res.resname] += 1
            continue
        if lastres is not None and lastres.resname not in bonddict:
            lastres = None
        # sometimes residues are next to each other in the chain,
        # but not in reality
        # we skip those
        if lastres is not None and res_dist(res, lastres) > 1:
            #~ print('lastres too far:', lastres, res, res_dist(res, lastres))
            #~ try:
                #~ print("Distance:", res['N'] - lastres['C'])
            #~ except KeyError:
                #~ pass
            lastres = None
            
        if nextres is not None and nextres.resname not in bonddict:
            #~ print('not found', nextres.resname)
            nextres = None
        if nextres is not None and res_dist(nextres,res) > 1:
            print('nextres too far:', res, nextres, res_dist(nextres,res))
            try:
                print("Distance:", nextres['N'] - res['C'])
            except KeyError:
                pass
            nextres = None
            
        cur_bonds, cur_backbone, err = getlengths(res, lastres)
        cur_angles, cur_backangles, angerr = getangles(res, lastres, nextres)
        if permissive or not err:
            this_bonds = res_bonds[res.resname]
            for k,v in list(cur_bonds.items()):
                this_bonds[k].append(v)
            for k,v in list(cur_backbone.items()):
                backbone_bonds[k].append(v)
        if permissive or not angerr:
            this_angles = res_angles[res.resname]
            for k,v in list(cur_angles.items()):
                this_angles[k].append(v)
            for k,v in list(cur_backangles.items()):
                back_angles[k].append(v)
        #~ print(lastres, res, nextres, len(back_angles[('CA','C','N')]), ('CA','C','N') in cur_backangles)


def applytobonds(func, *args, **kwargs):
    newres = dict()
    newangles = dict()
    for res, curdict in list(res_bonds.items()):
        newres[res] = {k:func(v, *args, **kwargs) for k,v in list(curdict.items())}
    newbdict = {k:func(v, *args, **kwargs) for k,v in list(backbone_bonds.items())}
    for res, curdict in list(res_angles.items()):
        newangles[res] = {k:func(v, *args, **kwargs) for k,v in list(curdict.items())}
    newbackang = {k:func(v, *args, **kwargs) for k,v in list(back_angles.items())}
    
    return newres, newbdict, newangles, newbackang

def stats(l):
    import numpy as np
    d = dict()
    d['mean'] = mean = np.mean(l)
    d['max'] = max(l)
    d['min'] = min(l)
    d['std'] = std = np.std(l)
    d['range'] = rnge = d['max'] - d['min']
    d['range/std'] = rnge / std
    d['range/mean'] = rnge / mean
    d['std/mean'] = std / mean
    d['N'] = len(l)
    return d

def makehist(l, nbins):
    from numpy import histogram
    return histogram(l,nbins)
    
def singleplot(res, hist, bins, *args, **kwargs):
    from matplotlib import pyplot
    xvals = (bins[1:] + bins[:-1])/2
    pyplot.plot(xvals, hist, *args, **kwargs)

def load_hists():
    import pickle as pickle
    return pickle.load(histfile.open('rb'))

def load_stats():
    import pickle as pickle
    return pickle.load(statfile.open('rb'))

if __name__ == '__main__':
    mydir = Dir('/home/wendell/idp/allpdbs')
    fnames = [p for p in mydir.children() if p.extension == 'ent']

    for n, fn in enumerate(fnames):
        print(n, fn)
        parsepdb(str(fn), True)
        #~ if n >= 10:
            #~ break

    print('Calculating statistics...')

    resstats, backstats, angstats, backangstats = applytobonds(stats)
    bondstats = ([(n, d) for n, d in list(resstats.items())]  + [('backbone',backstats)])
    
    res_hists, back_hists, ang_hists, back_ang_hists = applytobonds(makehist, 200)

    for name,d in bondstats:
        for bond, stat in list(d.items()):
            bname = bond[0] + '-' + bond[1]
            if stat['std/mean'] > .02:
                print(name, bname, 'std/mean', stat['std/mean'])
            if stat['range/mean'] > 1:
                print(name, bname, 'range/mean', stat['range/mean'])
    
    if True:
        print('Writing to file...')
        import pickle as pickle
        with statfile.open('wb') as f:
            pickle.dump((resstats, backstats, angstats, backangstats), f, -1)
            
        with histfile.open('wb') as f:
            pickle.dump((res_hists, back_hists, ang_hists, back_ang_hists), f, -1)
        print('finished.')
        print('Strange residues:', ", ".join([k+':'+str(v) for k,v in list(notfound.items())]))
        print(len(broken), 'Missing bonds.')
