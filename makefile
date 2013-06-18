UNAME := $(shell uname)
#CXX=${CXX}
SWIG=swig -Wextra -shadow -python -py3 -c++
CCOPTS=-Wall -O2 -fPIC

#INC=-I/usr/include/python2.7
INC=`python3-config --includes`

#~ ifeq ("$(UNAME)","Darwin")
	#~ CC=/opt/local/bin/g++-mp-4.6
	#~ SWIG=/opt/local/bin/swig
	#~ INC=-I/opt/local/Library/Frameworks/Python.framework/Versions/Current/include/python2.7
	#~ LIB=/opt/local/lib/libpython2.7.dylib
#~ endif

.PHONY: all 2d 3d printout clean

all: printout 2d 3d
	@echo "making all."

2d: _sim2d.so

3d: _sim.so

	
printout:
	@echo Running on \"$(UNAME)\"

clean:
	rm -f *.o *.so *.gch sim_wrap*.cxx

sim_wrap2d.cxx: sim.i collection.hpp constraints.hpp interaction.hpp vecrand.hpp collection.cpp constraints.cpp interaction.cpp vecrand.cpp vec.hpp
	$(SWIG) -DVEC2D sim.i
	mv sim_wrap.cxx sim_wrap2d.cxx

sim_wrap3d.cxx: sim.i collection.hpp constraints.hpp interaction.hpp vecrand.hpp collection.cpp constraints.cpp interaction.cpp vecrand.cpp vec.hpp
	$(SWIG) -DVEC3D sim.i
	mv sim_wrap.cxx sim_wrap3d.cxx

sim_wrap3d.o: sim_wrap3d.cxx
	$(CXX) $(CCOPTS) -DVEC3D -c sim_wrap3d.cxx $(INC)

sim_wrap2d.o: sim_wrap2d.cxx
	$(CXX) $(CCOPTS) -DVEC3D -c sim_wrap2d.cxx $(INC)

_sim2d.so: sim_wrap2d.o
	$(CXX) $(CCOPTS) -shared sim_wrap2d.o -o _sim2d.so $(LIB)
	
_sim.so: sim_wrap3d.o
	$(CXX) $(CCOPTS) -shared sim_wrap3d.o -o _sim.so $(LIB)
