import sys
import argparse
from array import array as pyarray

import numpy as np

from .statistics import StatSet
from . import util

class Simulation:
    def __init__(self, collec, dt, time, printn=100, cut=0.0):
        """
        cut is expected to be a fraction (e.g., 0.5), or (if larger than 1) assumed to be
        the amount of time to cut.
        """
        ndim = len(box.boxshape())
        if ndim == 2:
            import pyparm.d2 as sim
            self.sim = sim
        elif ndim == 3:
            import pyparm.d3 as sim
            self.sim = sim
        else:
            raise NotImplementedError("Unknown number of dimensions")

        self.collec = collec
        self.statsets = []

        self.dt = dt
        self.steps_done = 0 # steps completed
        self.steps_total = int(np.round(time / dt)) # total steps

        self.print_tot = printn
        self.printn = 0
        self.cut = cut

        self._progress = None

    @property
    def progress(self):
        if self._progress is None:
            self._progress = util.Progress(self.steps_total)
        return self._progress

    @property
    def cut(self):
        return self._cut
    
    @cut.setter
    def cut(self, value):
        if value <= 0:
            value = 0
        elif value <= 1: 
            value = int(np.round(value * self.steps_total)) # total steps
        else:
            value = int(np.round(value))
        self._cut = value

    def add_interaction(self, inter):
        self.collec.addInteraction(inter)
    
    def add_tracker(self, tracker):
        self.collec.addTracker(tracker)
    
    def add_stats(self, statset, statdt=None):
        if statdt is not None:
            statset.statdt = statdt 
        else:
            assert statset.statdt != None
        self.statsets.append(statset)

    def output(self, *args, **kwargs):
        print(*args, **kwargs)
        f = kwargs.get('file', sys.stdout)
        f.flush()

    def equilibrate(self, progress=True):
        prog = self.progress # this initializes self._progress
        for t in range(self.steps_done, self.cut):
            if t > 0: collec.timestep()
            self.steps_done = t
            if progress: self.progress_out()
            

    def progress_out(self, force=False):
        if self.print_tot <= 0: return
        t = self.steps_done
        printt = float(self.steps_total) * self.printn / self.print_tot
        if force or t >= printt:
            self.output('{9:.6g} ---- '.format(t * self.dt), prog.eta_str(t))
            while t >= printt:
                self.printn += 1
                printt = float(self.steps_total) * self.printn / self.print_tot

    def progress_str(self, time):
        return '{9:.6g} ---- '.format(time)

    def run(self, progress=True, print_updates=False):
        self.equilibrate(progress=progress)
        statts = [self.steps_done for _ in self.statsets]
        
        for t in range(self.steps_done, self.steps_total):
            if t > 0: collec.timestep()
            self.steps_done = t

            for n, stime in enumerate(statts):
                if t >= stime:
                    statset = self.statsets[n]
                    if statset.statdt <= 0:
                        # ignore these
                        statts[n] = stime = self.steps_total
                        continue
                    while t >= stime:
                        stime += statset.statdt
                        statts[n] = stime
                    if print_updates: print('Updating', statset)
                    statset.update(t * self.dt)

            if progress: self.progress_out()

        for _, statset in self.statsets:
            statset.safe_write()

    @staticmethod
    def arg_parser(parser=None, **kw):
        if parser is None:
            kw.setdefault('formatter_class', argparse.ArgumentDefaultsHelpFormatter)
            parser = argparse.ArgumentParser(**kw)
        group = parser.add_argument_group('Global Simulation Settings')

        group.add_argument('-d', '--dt', type=float, default=0.2, help='timestep')
        group.add_argument('-t', '--time', type=int, default=2000, help='total run time')
        group.add_argument('--printn', type=int, default=200, help='Number of times to print progress')

        group.add_argument('-x', '--cut', type=float, default=0.5, help='Amount of data to cut')

        group.add_argument('-O', '--outfilename', default='test/t{time}', help='base path name for output files')
        return parser