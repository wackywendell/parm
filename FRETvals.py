from numpy import sqrt, mean, std, array

Forsterdist = 54

ijset = [(9,130), (33,130), (54,130), (72,130), (33,72), (9,54),
        (72,92), (54,72), (9,72), (9,33), (54,92), (92,130)]

ETeffs = {
          (9,130) : 0.36, (33,130) : 0.53,(54,130) : 0.5, (72,130) : 0.55, (92,130) : 0.65, (33,72) : 0.72,
          (9,54) : 0.7, (72,92) : 0.85, (54,72) : 0.85, (9,72) : 0.64, (9,33) : 0.86, (54,92) : 0.66}
          
ETerrs = { (9,130) : 0.0529, (33,130) : 0.0608, (54,130) : 0.0627,
            (72,130) : 0.0621, (92,130) : 0.0624, (33,72) : 0.0643,
            (9,54) : 0.0608, (72,92) : 0.0451, (54,72) : 0.0437,
            (9,72) : 0.0658, (9,33) : 0.0503, (54,92) : 0.0674 }

ijs = [(54, 72), (72, 92), (9, 33), (54, 92), (92, 130), (33, 72),
            (9, 54), (72, 130), (9, 72), (54, 130), (33, 130), (9, 130)]
ETs = [ETeffs[ij] for ij in ijs]
ETerrlst = [ETerrs[ij] for ij in ijs]

ETeffs3 = {
          (9,130) : 0.36, (33,130) : 0.70,(54,130) : 0.67, (72,130) : 0.77, (92,130) : 0.86, (33,72) : 0.67,
          (9,54) : 0.70, (72,92) : 0.93, (54,72) : 0.91}
ijs3 = [ij for ij in ijs if ij in ETeffs3]
ETs3 = [ETeffs3[ij] for ij in ijs3]
ETerrlst3 = [ETerrs[ij] for ij in ijs3]

revijs = sorted(ijset, key=lambda i_j: -abs(i_j[1]-i_j[0]))
expETsij = [ETeffs[ij] for ij in revijs]
expETerrs = [ETerrs[ij] for ij in revijs]

ijexp = sorted(ijset, key=lambda ij: -ETeffs[ij])
expETs = [ETeffs[ij] for ij in ijexp]

ijsRW = [(72, 92), (54, 72), (9, 33), (54, 92), (92, 130), (33, 72),
            (9, 54), (72, 130), (9, 72), (54, 130), (33, 130), (9, 130)]
expETsRW = [ETeffs[ij] for ij in ijsRW]
expETerrsRW = [ETerrs[ij] for ij in ijsRW]

FRETdists={
         (9,130) : (63.4,16.3), (33,130) : (52.8,13.2), (54,130) : (54.0,13.5), (72,130) : (51.6,12.9),
         (92,130) : (47.0,11.5), (33,72) : (44.0,10.5), (9,54) : (44.6,10.7), (72,92) : (37.9,8.5),
         (54,72) : (37.9,8.5), (9,72) : (47.7,11.5), (9,33) : (37.3,8.3), (54,92) : (46.8,11.1)}

#def ETeffrsq(ds):
#    dists = [(d - ed)**2 for d,ed in zip(ds, expETs)]
#    return sqrt(mean(dists))

def ETeffrsq(d, ph3=False):
    ETdict = ETeffs if not ph3 else ETeffs3
    ijset = ijs if not ph3 else ijs3
    dists = [(d[ij] - ETdict[ij])**2 for ij in ijset]
    return sqrt(mean(dists))

def ETeffrsqerrs(simETs):
    dxs = array([(simETs[ij] - ETeffs[ij]) for ij in ijs])
    sigys = array([ETerrs[ij] for ij in ijs])
    N = len(ETeffs)
    
    D = sqrt(mean(dxs**2))
    
    dDdxi = dxs/D/N
    Derr = sqrt(sum((dDdxi * sigys)**2))
    
    return D, Derr

def __ETeffrsqerrs_old(simETs):
    dists = [(simETs[ij] - ETeffs[ij])**2 for ij in ETeffs]
    D = sqrt(mean(dists))
    disterrs = [(ETerrs[ij] * abs(simETs[ij] - ETeffs[ij])/D/12.0)**2
                        for ij in ETerrs]
    Derr = sqrt(sum(disterrs))
    return D, Derr

distribijs = [(1, 140), (12, 127), (14, 125), (21, 118), (21, 111), (36, 75), (38, 77), (43, 82), (57, 133), (59, 135), (63, 139), (9, 130), (33, 130), (54, 130), (72, 130), (92, 130), (33, 72), (9, 54), (72, 92), (54, 72), (9, 72), (9, 33), (54, 92)]

SEQ = "\
MET ASP VAL PHE MET LYS GLY LEU SER LYS ALA LYS GLU \
GLY VAL VAL ALA ALA ALA GLU LYS THR LYS GLN GLY VAL \
ALA GLU ALA ALA GLY LYS THR LYS GLU GLY VAL LEU TYR \
VAL GLY SER LYS THR LYS GLU GLY VAL VAL HIS GLY VAL \
ALA THR VAL ALA GLU LYS THR LYS GLU GLN VAL THR ASN \
VAL GLY GLY ALA VAL VAL THR GLY VAL THR ALA VAL ALA \
GLN LYS THR VAL GLU GLY ALA GLY SER ILE ALA ALA ALA \
THR GLY PHE VAL LYS LYS ASP GLN LEU GLY LYS ASN GLU \
GLU GLY ALA PRO GLN GLU GLY ILE LEU GLU ASP MET PRO \
VAL ASP PRO ASP ASN GLU ALA TYR GLU MET PRO SER GLU \
GLU GLY TYR GLN ASP TYR GLU PRO GLU ALA".split(' ')
