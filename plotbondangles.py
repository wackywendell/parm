from matplotlib import pylab as plt
import lengthanalysis as lengths

bonds, bbonds, angles, bangles = lengths.load_hists()
bondstat, bbondstat, anglestat, banglestat = lengths.load_stats()

def group_plot(hists, stats=None, title=None, show=True, lims=True):
    xmin = None
    xmax = None
    for k,v in sorted(hists.items()):
        label="--".join(k)
        if stats:
            try:
                mean = stats[k]['mean']
            except KeyError:
                print((title, k, list(stats.keys())))
                raise
            std = stats[k]['std']
            if xmax is None or xmax < mean + 3*std:
                xmax =  mean + 3*std
            if  xmin is None or xmin > mean - 3*std:
                xmin =  mean - 3*std
            label += ", $\mu=%.2f,$ $\sigma=%.3f$ $(%.1f\%%)$" % (mean, std, std/mean * 100)
            
        lengths.singleplot(k,*v, label=label)
        
    if len(hists) > 0: plt.legend(loc=0)
    #~ print(xmin, xmax)
    if xmin and xmax and lims:
        plt.xlim(xmin,xmax + (xmax-xmin)*.0001)
    else:
        plt.xlim(0,None)
    title and plt.title(title)
    show and plt.show()


def plot_each(hists, stats=None, savefname = None, show = False, *args, **kwargs):
    plt.figure(figsize=(15,12))
    for k,v in list(hists.items()):
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
        show and plt.show()

allbonds = dict(bonds)
allbonds['backbone'] = bbonds
allbondstat = dict(bondstat)
allbondstat['backbone'] = bbondstat
allangles = dict(angles)
allangles['backbone'] = bangles
allanglestat = dict(anglestat)
allanglestat['backbone'] = banglestat

#~ group_plot(bangles, banglestat)
print('bonds 1')
plot_each(allbonds, allbondstat, '/home/wendell/idp/blengths/%s-bonds.eps')
print('bonds 2')
plot_each(allbonds, allbondstat, '/home/wendell/idp/blengths/%s-bonds-0.eps', lims=False)

print('angles 1')
plot_each(allangles, allanglestat, '/home/wendell/idp/blengths/%s-angles.eps')
print('angles 2')
plot_each(allangles, allanglestat, '/home/wendell/idp/blengths/%s-angles-0.eps', lims=False)
