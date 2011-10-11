CC = gcc
CC += -c
CPP = g++ 
CPP += -c -Wall -O2 -fPIC
LD = g++ -Wall -O2
SO = g++ -Wall -shared -O2

b = ../bin
ALL = interaction.o collection.o _sim.so
#simulation.so

all: $(ALL)

PYOPT = -O2 -fPIC -I/usr/include/python2.7 -I/usr/include
PYSOOPT = -O2 -shared -Wl,--export-dynamic -Wl,-no-undefined -L/usr/lib -lboost_python -L/usr/lib/python2.7/config -lpython2.7

#FULL = $b/test1 ...
#full: $(FULL)

clean:
	rm -f *.o *.so $(ALL) *.gch sim_wrap.cxx

interaction.o: interaction.hpp interaction.cpp vec.hpp
	$(CPP) $^

collection.o: vec.hpp interaction.hpp collection.hpp collection.cpp vecrand.cpp vecrand.hpp
	$(CPP) $^

simulation.o: vec.hpp interaction.hpp collection.hpp simulation.cpp
	g++ $(PYOPT) -c simulation.cpp

simulation.so: simulation.o collection.o interaction.o
	g++ $(PYSOOPT) -o $@ $^

sim_wrap.cxx: vec.hpp interaction.hpp collection.hpp interaction.cpp collection.cpp vecrand.cpp vecrand.hpp sim.i
	swig -shadow -python -c++ sim.i

sim_wrap.o: sim_wrap.cxx
	g++ -fPIC -c sim_wrap.cxx -I/usr/include/python2.7

_sim.so: sim_wrap.o
	g++ -shared sim_wrap.o -o _sim.so

#$b/vectest: vectest.cpp vec.hpp
#	$(LD) -o $@ $^

%.o: %.cpp %.hpp
	$(CPP) $^

#molecule.so: molecule.o heptanemol.o
#	$(SO) -o $@ $^
