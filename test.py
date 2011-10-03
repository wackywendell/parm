import sim

a = sim.atomvec(5, sim.fvector([3,1,3]))

i = sim.spring(3, 1)

ip = sim.intraMolNNPair(sim.avector([a]),i)
