import xyzstats
from xyzstats import statkeeper

smFRETs={
(9,130) : (63.4,16.3),
(33,130) : (52.8,13.2),
(54,130) : (54.0,13.5),
(72,130) : (51.6,12.9),
(92,130) : (47.0,11.5),
(33,72) : (44.0,10.5),
(9,54) : (44.6,10.7),
(72,92) : (37.9,8.5),
(54,72) : (37.9,8.5),
(9,72) : (47.7,11.5),
(9,33) : (37.3,8.3),
(54,92) : (46.8,11.1)
}

def plot_reslengths(sk):
    Rijs = [np.array(sk.Rij(i,j), dtype=float) for i,j in sorted(smFRETs.keys())]
    vals = [(Rij.mean(), Rij.std()) for Rij in Rijs]
    means, stds = list(zip(*vals))
    
        
