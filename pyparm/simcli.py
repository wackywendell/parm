import sys
import argparse
from array import array as pyarray

import numpy as np

from .statistics import StatSet
from . import util

class Simulation:
    def __init__(self, atoms, box, collec, dt, time, printsteps=1000, cut=0.0):
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

        self.atoms = atoms
        self.box = box
        self.collec = collec
        self.interactions = []
        self.trackers = []
        self.statsets = []

        self.dt = dt
        self.steps_done = 0 # steps completed
        self.steps_total = int(np.round(time / dt)) # total steps

        self.printsteps = printsteps
        self.printn = 0
        if cut <= 0:
            self.cut = 0
        elif cut <= 1: 
            self.cut = int(np.round(cut * self.steps_total)) # total steps
        else:
            self.cut = int(np.round(cut))

        self._progress = None

    @property
    def progress(self):
        if self._progress is None:
            self._progress = util.Progress(self.steps_total)
        return self._progress

    def add_interaction(self, inter):
        self.collec.addInteraction(inter)
        self.interactions.append(inter)
    
    def add_tracker(self, tracker):
        self.collec.addTracker(tracker)
        self.trackers.append(tracker)

    def add_stats(self, statset, statdt=None):
        if statdt is None:
            statdt = statset.statdt
        self.statsets.append((statdt, statset))

    def output(self, *args, **kwargs):
        print(*args, **kwargs)
        f = kwargs.get('file', sys.stdout)
        f.flush()

    def equilibrate(self, progress=True):
        prog = self.progress # this initializes self._progress
        printt = float(self.steps_total) * self.printn / self.printsteps
        for t in range(self.steps_done, self.cut):
            if t > 0: collec.timestep()
            self.steps_done = t
            if progress: self.progress_out()

            

    def progress_out(self, force=False):
        t = self.steps_done
        printt = float(self.steps_total) * self.printn / self.printsteps
        if force or t >= printt:
            self.output('{9:.6g} ---- '.format(t * self.dt), prog.eta_str(t))
            while t >= printt:
                self.printn += 1
                printt = float(self.steps_total) * self.printn / self.printsteps

    def progress_str(self, time):
        return '{9:.6g} ---- '.format(time)

    def run(self, progress=True):
        self.equilibrate(progress=progress)
        printt = float(self.steps_total) * self.printn / self.printsteps
        statts = [self.steps_done for _ in self.statsets]
        
        for t in range(self.steps_done, self.steps_total):
            if t > 0: collec.timestep()
            self.steps_done = t

            for n, stime in enumerate(statts):
                if t >= stime:
                    sdt, statset = self.statsets[n]
                    while t >= stime:
                        stime += sdt
                        statts[n] = stime
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

