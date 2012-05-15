UNAME := $(shell uname)
#CXX=${CXX}
SWIG=swig
CCOPTS=-Wall -O3 -fPIC

#INC=-I/usr/include/python2.7
INC=`python2-config --includes`

#~ ifeq ("$(UNAME)","Darwin")
	#~ CC=/opt/local/bin/g++-mp-4.6
	#~ SWIG=/opt/local/bin/swig
	#~ INC=-I/opt/local/Library/Frameworks/Python.framework/Versions/Current/include/python2.7
	#~ LIB=/opt/local/lib/libpython2.7.dylib
#~ endif

ALL=_sim.so
all: printout $(ALL)

printout:
	@echo Running on \"$(UNAME)\"

clean:
	rm -f *.o *.so $(ALL) *.gch sim_wrap.cxx

sim_wrap.cxx: *.hpp *.cpp sim.i
	$(SWIG) -Wextra -shadow -python -c++ sim.i

sim_wrap.o: sim_wrap.cxx
	$(CXX) $(CCOPTS) -c sim_wrap.cxx $(INC)

_sim.so: sim_wrap.o
	$(CXX) $(CCOPTS) -shared sim_wrap.o -o _sim.so $(LIB)
