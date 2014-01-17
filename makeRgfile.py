from optparse import OptionParser
import xyzfile
import numpy as np
import FRETs

parser = OptionParser()
parser.add_option('-n', type=int, default=0)
parser.add_option('-s', dest='skip', type=int, default=1)
parser.add_option('-p', type=int, default=5)

Rijkeys = [k[:-3] for k in dir(FRETs) if k.endswith('ijs')] + ['auto', 'none']
parser.add_option('-R', '--Rijs', type='choice', choices=Rijkeys, default='auto')

opts, args = parser.parse_args()
fname = args[0]
bname, _, _ = fname.rpartition('.')
bbname = fname.split('/')[-1]

if opts.Rijs == 'auto':
    prot = bbname.split('-')[0].lower()
    lowkeys = {k.lower():k for k in Rijkeys}
    if prot not in lowkeys:
        print('Rij key not found', repr(prot))
        opts.Rijs = None
    else:
        opts.Rijs = lowkeys[prot]

if opts.Rijs is not None and opts.Rijs.lower() != 'none':
    ijs = getattr(FRETs, opts.Rijs + 'ijs')
else:
    ijs = []

#~ Rgf = bname + '-Rgs.tsv'
#~ print('saving to', Rgf)
Rgz = bname + '.rgs.npz'
print('saving to', Rgz)
if ijs:
    rijz = bname + '.rijs.npz'
    print('saving to', rijz)

def Rg_from_arr(arr):
    com = np.mean(arr, 0)
    return np.sqrt(np.sum(np.mean((arr-com)**2, 0)))

def get_Rg_ns(locs):
    global arr
    arr = locs
    ns, Rgs = np.arange(2, len(locs)+1), []
    for n in ns:
        Rgs.append([Rg_from_arr(locs[m:n+m]) for m in range(0, len(locs)-int(n)+1)])
    return ns, Rgs
    

allns, allRgs, allRgstd = [], [], []
allrijs, allts = [], []
def write_npzs():
    global allns, allRgs, allRgstd
    global ijs, allrijs, allts
    #~ np.savetxt(str(Rgf), np.array([curns, curRgs]).T, fmt='%.18g', delimiter='\t')
    curns = np.mean(allns, 0)
    curRgs = np.mean(allRgs, 0)
    Rgstd = np.std(allRgs, 0)

    np.savez_compressed(Rgz, ns=curns, Rgs=curRgs, Rgstd=Rgstd)
    if not ijs: return
    rijs = np.array(allrijs)
    ETeffs = 1.0 / (1.0 + (rijs**6.0 / FRETs.Forsterdist**6.0))
    np.savez_compressed(rijz, ijs=ijs, rijs=rijs, ET=np.mean(ETeffs,0), ETstd = np.std(ETeffs,0), ts=allts)
    
with open(fname) as f:
    xyz = xyzfile.XYZreader(f)
    frames = iter(xyz)
    for fr in frames:
        ns, Rgs = get_Rg_ns(fr.locarray)
        allns.append(np.array(ns))
        allRgs.append(np.array([np.mean(l) for l in Rgs]))
        allRgstd.append(np.array([np.std(l) for l in Rgs]))
                
        if ijs:
            currijs = np.sqrt(np.array(
                [np.sum((fr.locarray[i-1] - fr.locarray[j-1])**2) for i,j in ijs]
            ))
            allrijs.append(currijs)
            allts.append(fr.time)
            
        
        if opts.p > 0 and len(allns) % opts.p == 0:
            print(len(allns), int(fr.time))
            write_npzs()
        
        if opts.n > 0 and len(allns) >= opts.n:
            break
        for _ in range(opts.skip - 1):
            #print('skip', fr.time)
            try:
                fr = next(frames)
            except StopIteration:
                break
            continue

write_npzs()
#~ ns = np.mean(allns, 0)
#~ Rgs = np.mean(allRgs, 0)
#~ 
#~ np.savetxt(str(Rgf), np.array([ns, Rgs]).T, fmt='%.18g', delimiter='\t')
#~ np.savez_compressed(Rgz, ns=ns, Rgs=allRgs, Rgstd=allRgstd)
