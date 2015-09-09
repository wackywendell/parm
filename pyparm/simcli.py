# -*- coding: utf-8 -*-

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
        self.collec.add_interaction(inter)
    
    def add_tracker(self, tracker):
        self.collec.add_tracker(tracker)
    
    def add_stats(self, statset):
        assert statset.statdt != None
        self.statsets.append(statset)

    @property
    def time(self):
        return self.dt * self.steps_done
    
    def output(self, *args, **kwargs):
        print(*args, **kwargs)
        f = kwargs.get('file', sys.stdout)
        f.flush()

    def equilibrate(self, progress=True):
        prog = self.progress # this initializes self._progress
        for t in range(self.steps_done, self.cut+1):
            if t > 0: self.collec.timestep()
            self.steps_done = t
            if progress: self.progress_out()
            

    def progress_out(self, force=False):
        if self.print_tot <= 0: return
        t = self.steps_done
        printt = float(self.steps_total) * self.printn / self.print_tot
        if force or t >= printt:
            self.output(self.progress_str(t * self.dt), self.progress.eta_str(t))
            while t >= printt:
                self.printn += 1
                printt = float(self.steps_total) * self.printn / self.print_tot

    def progress_str(self, time):
        return '{:9.6g} ---- '.format(time)

    def run(self, progress=True, print_updates=False):
        self.equilibrate(progress=progress)
        # t, statts, stime are all in units of steps
        statts = [self.steps_done for _ in self.statsets]
        
        startt = self.steps_done
        for t in range(self.steps_done, self.steps_total+1):
            if t > startt: self.collec.timestep()
            self.steps_done = t

            for n, stime in enumerate(statts):
                if t >= stime:
                    statset = self.statsets[n]
                    if statset.statdt <= 0:
                        # ignore these
                        statts[n] = stime = self.steps_total
                        continue
                    while t >= stime:
                        stime += statset.statdt / self.dt
                        statts[n] = stime
                    if print_updates: self.output('Updating', statset, 'time', t, 'dt', statset.statdt)
                    statset.update(t * self.dt)

            if progress: self.progress_out()
            
        for statset in self.statsets:
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
