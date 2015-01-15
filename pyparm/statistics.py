import argparse
from collections import OrderedDict
from array import array as pyarray

import numpy as np

class Statistic:
    """a class that updates itself with gather(), and will return
    a numpy Array (or dict of numpy arrays) with stats()"""
    def gather(self, time):
        raise NotImplementedError

    def stats(self):
        raise NotImplementedError

def get_order(atoms, box, local=True, weighted=True):
    import tess
    locs = [tuple(a.x) for a in atoms]
    limits = tuple(box.boxshape())
    cntr = tess.Container(locs, limits, radii=sigmas/2., periodic=True)
    return cntr.order(local=local, weighted=weighted)

class SimpleStat(Statistic):
    """A statistic that is a simple function of time. This class calls a given function
    at each `gather()` time, and manages the array of previous values."""
    def __init__(self, func, name=None):
        """
        func: a function that returns a simple Statistic when run as func(atoms, box, collection, time)
        name: the name of the statistic; defaults to function name"""
        self.func = func
        self.arr = []
        self.name = func.__name__ if name is None else name

    def gather(self, time):
        values = self.func(time)
        if len(self.arr) == 0:
            if np.shape(values) != ():
                self.arr = [pyarray('d') for _ in values]
                self.arr_flat = False
            else:
                self.arr = pyarray('d')
                self.arr_flat = True
        
        if self.arr_flat:
            self.arr.append(values)
        else:
            for v,a in zip(values, self.arr):
                a.append(v)

    def stats(self):
        return np.array(self.arr, dtype=float)

    @classmethod
    def get_Ql(cls, atoms, box, weighted=True, name=None):
        def Ql(time):
            return get_order(atoms, box, local=True, weighted=weighted)
        return cls(Ql, name=name)
    
    @classmethod
    def get_Qg(cls, atoms, box, weighted=True, name=None):
        def Qg(time):
            return get_order(atoms, box, local=False, weighted=weighted)
        return cls(Qg, name=name)

    @classmethod
    def get_events(cls, collec, name=None):
        def events(time):
            return collec.events_processed
        return cls(events, name=name)

class ReturnStat(Statistic):
    """A statistic that is self-updating, and only needs to call a simple function to get its values.

    Any `statetracker`, for example, probably falls into this category."""
    def __init__(self, func, name=None):
        """
        func: a function that returns a simple Statistic when run as func()
        name: the name of the statistic; defaults to function name"""
        self.func = func
        self.name = func.__name__ if name is None else name

    def gather(self, time):
        """This statistic does not need to be told to gather; it does it on its own."""
        pass

    def stats(self):
        return self.func()

class FixedStat(Statistic):
    """
    A class that simply returns the same thing every time.
    """
    def __init__(self, val, name):
        self.val = np.array(val)
        self.name = name

    def gather(self, time):
        """This statistic does not change."""
        pass

    def stats(self):
        return self.val

class StatSet(OrderedDict):
    def __init__(self, fname, statdt, writestep=20):
        OrderedDict.__init__(self)
        self.fname = fname
        self.writen = 0
        self.writestep = writestep
        self.print_on_write = False
        self.statdt = statdt
        self.writes = 0

    def add(self, stat):
        assert stat.name not in self
        self[stat.name] = stat

    def add_func(self, func, name=None):
        s = SimpleStat(func, name=name)
        self.add(s)
        return func

    def add_return(self, func, name=None):
        s = ReturnStat(func, name=name)
        self.add(s)
        return func

    def add_fixed(self, dtype=None, **kw):
        stats = {}
        for k,v in kw.items():
            stat = stats[k] = FixedStat(np.array(v, dtype=dtype), name=k)
            self.add(stat)
        return stats

    def update(self, time, write=None):
        for name,s in self.items():
            s.gather(time)

        self.writen += 1
        self.writen %= max(self.writestep, 1)
        if write or (write is None and self.writen == 0):
            print('Writing at', time)
            self.safe_write()

    def add_basics(self, atoms, box, collec):
        """Add the basic stats: [t, P, U, K, E, T, COMV, V]"""
        @self.add_func
        def t(time):
            return time

        @self.add_func
        def P(time):
            return collec.pressure()

        @self.add_func
        def U(time):
            return collec.potentialenergy()

        @self.add_func
        def K(time):
            return collec.kinetic()

        @self.add_func
        def E(time):
            return collec.energy()

        @self.add_func
        def T(time):
            return collec.temp()

        @self.add_func
        def COMV(time):
            return collec.comv().mag()

        @self.add_func
        def V(time):
            return box.V()

    def get_statistics(self):
        return {name:s.stats() for name,s in self.items()}

    def safe_write(self):
        if self.print_on_write:
            print('Writing to', self.fname)
        try:
            np.savez_compressed(str(self.fname), **self.get_statistics())
            self.writes += 1
        except:
            print("Interrupted in save, retrying...")
            np.savez_compressed(str(self.fname), **self.get_statistics())
            self.writes += 1
            raise

    @staticmethod
    def arg_parser(*stats, parser=None, title=None, description=None, defaultn=1000):
        if parser is None:
            parser = argparse.ArgumentParser(title=title, description=description)
        for stat in stats:
            group = parser.add_argument_group('{} Settings'.format(stat))

            tgroup = group.add_mutually_exclusive_group()
            tgroup.add_argument('--{}time'.format(stat), type=float, help='amount of time between taking statistics')
            tgroup.add_argument('--{}n'.format(stat), type=int, default=defaultn, help='Number of statistics to take')

            group.add_argument('--{}writen'.format(stat), type=int, default=20, help='Number of updates between disk writes')

        def new_parse_args(*args, **kwargs):
            parsed_args = argparse.ArgumentParser.parse_args(parser, *args, **kwargs)
            if 'statistic_sets' not in parsed_args:
                parsed_args.statistic_sets = stats
            return parsed_args

        parser.parse_args = new_parse_args

    @classmethod
    def from_opts(cls, opts):
        """
        """
        statsets = []
        stats = opts.statistic_sets
        
        for stat in stats:
            fname = opts.outfilename + '.' + str(stat) + '.npz'
            fname = fname.format(**vars(opts))

            statdt = vars(opts).get('{}time'.format(stat))
            if statdt is None:
                statn = vars(opts).get('{}n'.format(stat))
                statdt = (1 - opts.cut) * opts.time / statn

            statset = cls(fname=fname, writestep=vars(opts).get('{}writen'.format(stat)))
            statset.statdt = statdt
            statsets.append(statset)
        return statsets

    def __repr__(self):
        return "StatSet({})".format(repr(self.fname))

    def __str__(self):
        lst = [self.fname] + [str(k) for k in self]
        return "StatSet(" + ', '.join(lst) + ')'