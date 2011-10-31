from matplotlib import pylab as plt
import lengthanalysis as lengths

bonds, bbonds, angles, bangles = lengths.load_hists()
bondstat, bbondstat, anglestat, banglestat = lengths.load_stats()

def group_plot(hists, stats=None, title=None, show=True, lims=True):
    xmin = None
    xmax = None
    for k,v in sorted(hists.items()):
        label="-".join(k)
        if stats:
            try:
                mean = stats[k]['mean']
            except KeyError:
                print(title, k, stats.keys())
                raise
            std = stats[k]['std']
            if xmax is None or xmax < mean + 3*std:
                xmax =  mean + 3*std
            if  xmin is None or xmin > mean - 3*std:
                xmin =  mean - 3*std
            label += ", $\sigma=%.2f\,(%.1f\%%)$" % (std, std/mean * 100)
            
        lengths.singleplot(k,*v, label=label)
        
    plt.legend(loc=0)
    #~ print(xmin, xmax)
    if xmin and xmax and lims:
        plt.xlim(xmin,xmax + (xmax-xmin)*.0001)
    else:
        plt.xlim(0,None)
    if title:
        plt.title(title)
    if show:
        plt.show()


def plot_each(hists, stats=None, savefname = None, *args, **kwargs):
    figure(figsize=(15,12))
    for k,v in hists.items():
        stat = None
        if stats:
            stat = stats[k]
        if savefname:
            plt.clf()
            #~ plt.gcf().set_size_inches((12,10),True)
        else:  plt.figure()
        group_plot(v, stat, *args, title=k, show=False, **kwargs)
        if savefname:
            plt.savefig(str(savefname) % k)
        plt.show()

allbonds = dict(bonds)
allbonds['backbone'] = bbonds
allbondstat = dict(bondstat)
allbondstat['backbone'] = bbondstat
allangles = dict(angles)
allangles['backbone'] = bangles
allanglestat = dict(anglestat)
allanglestat['backbone'] = banglestat

#~ group_plot(bangles, banglestat)
plot_each(allbonds, allbondstat, '/home/wendell/idp/blengths/%s-bonds.eps')
plot_each(allangles, allanglestat, '/home/wendell/idp/blengths/%s-angles.eps')

plot_each(allbonds, allbondstat, '/home/wendell/idp/blengths/%s-bonds-0.eps', lims=False)
plot_each(allangles, allanglestat, '/home/wendell/idp/blengths/%s-angles-0.eps', lims=False)
