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

ALL=_sim.so
2D=

all: printout $(ALL)
	@echo "making all."


2d: _sim2d.so
	
printout:
	@echo Running on \"$(UNAME)\"

clean:
	rm -f *.o *.so $(ALL) *.gch sim_wrap.cxx

sim_wrap.cxx: *.hpp *.cpp sim.i
	$(SWIG) $(2D) sim.i

sim_wrap.o: sim_wrap.cxx
	$(CXX) $(CCOPTS) -c sim_wrap.cxx $(INC)

2dln: vecrand2d.hpp vecrand2d.cpp
	[[ -e sim_wrap.cxx ]] && mv sim_wrap.cxx sim_wrap3d.cxx || true
	[[ -e sim_wrap.o ]] && mv sim_wrap.o sim_wrap3d.o || true
	[[ -e sim_wrap2d.cxx ]] && mv sim_wrap2d.cxx sim_wrap.cxx || true
	[[ -e sim_wrap2d.o ]] && mv sim_wrap2d.o sim_wrap.o || true
	ln -sf vecrand2d.hpp vecrand.hpp
	ln -sf vecrand2d.cpp vecrand.cpp
	$(eval 2D=-D2D)

_sim2d.so: 2dln sim_wrap.o
	$(CXX) $(CCOPTS) -shared sim_wrap.o -o _sim2d.so $(LIB)
	ln -sf vecrand3d.hpp vecrand.hpp
	ln -sf vecrand3d.cpp vecrand.cpp
	mv sim_wrap.cxx sim_wrap2d.cxx
	mv sim_wrap.o sim_wrap2d.o

_sim.so: sim_wrap.o
	$(CXX) $(CCOPTS) -shared sim_wrap.o -o _sim.so $(LIB)
