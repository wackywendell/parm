UNAME := $(shell uname)
#CXX=${CXX}
SWIG=swig -Wextra -shadow -python -py3 -c++
CCOPTS=-Wall -O2 -fPIC  -Wconversion

#INC=-I/usr/include/python2.7
INC=`python3-config --includes`

#~ ifeq ("$(UNAME)","Darwin")
	#~ CC=/opt/local/bin/g++-mp-4.6
	#~ SWIG=/opt/local/bin/swig
	#~ INC=-I/opt/local/Library/Frameworks/Python.framework/Versions/Current/include/python2.7
	#~ LIB=/opt/local/lib/libpython2.7.dylib
#~ endif

.PHONY: all 2d 3d 2dlong 3dlong printout clean wraps

all: 2d 2dlong 3d 3dlong
	@echo "making all."

2d: _sim2d.so

3d: _sim.so

2dlong: _sim2dlong.so

3dlong: _sim3dlong.so

wraps: sim_wrap2d.cxx sim_wrap2dlong.cxx sim_wrap3d.cxx sim_wrap3dlong.cxx

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

sim_wrap2d.o: sim_wrap2d.cxx
	$(CXX) $(CCOPTS) -DVEC2D -c sim_wrap2d.cxx $(INC)

sim_wrap3d.o: sim_wrap3d.cxx
	$(CXX) $(CCOPTS) -DVEC3D -c sim_wrap3d.cxx $(INC)

_sim2d.so: sim_wrap2d.o
	$(CXX) $(CCOPTS) -DVEC2D -shared sim_wrap2d.o -o _sim2d.so $(LIB)
	
_sim.so: sim_wrap3d.o
	$(CXX) $(CCOPTS) -DVEC3D -shared sim_wrap3d.o -o _sim.so $(LIB)

sim_wrap2dlong.cxx: sim.i collection.hpp constraints.hpp interaction.hpp vecrand.hpp collection.cpp constraints.cpp interaction.cpp vecrand.cpp vec.hpp
	$(SWIG) -DVEC2D -DLONGFLOAT sim.i
	mv sim_wrap.cxx sim_wrap2dlong.cxx

sim_wrap3dlong.cxx: sim.i collection.hpp constraints.hpp interaction.hpp vecrand.hpp collection.cpp constraints.cpp interaction.cpp vecrand.cpp vec.hpp
	$(SWIG) -DVEC3D -DLONGFLOAT sim.i
	mv sim_wrap.cxx sim_wrap3dlong.cxx

sim_wrap2dlong.o: sim_wrap2dlong.cxx
	$(CXX) $(CCOPTS) -DVEC2D -DLONGFLOAT -c sim_wrap2dlong.cxx $(INC)

sim_wrap3dlong.o: sim_wrap3dlong.cxx
	$(CXX) $(CCOPTS) -DVEC3D -DLONGFLOAT -c sim_wrap3dlong.cxx $(INC)

_sim2dlong.so: sim_wrap2dlong.o
	$(CXX) $(CCOPTS) -DVEC2D -DLONGFLOAT -shared sim_wrap2dlong.o -o _sim2dlong.so $(LIB)
	
_sim3dlong.so: sim_wrap3dlong.o
	$(CXX) $(CCOPTS) -DVEC3D -DLONGFLOAT -shared sim_wrap3dlong.o -o _sim3dlong.so $(LIB)

VECOPTS := 2D 3D

define VEC_TARGET_RULE
vecrand$(1).o: vec.hpp vecrand.hpp vecrand.cpp
	$(CXX) $(CCOPTS) -DVEC$(1) -c vecrand.cpp -o vecrand$(1).o

interaction$(1).o: vec.hpp vecrand.hpp interaction.hpp interaction.cpp
	$(CXX) $(CCOPTS) -DVEC$(1) -c interaction.cpp -o interaction$(1).o

constraints$(1).o: vec.hpp vecrand.hpp interaction.hpp constraints.hpp constraints.cpp
	$(CXX) $(CCOPTS) -DVEC$(1) -c constraints.cpp -o constraints$(1).o
		
collection$(1).o: vec.hpp vecrand.hpp interaction.hpp collection.hpp constraints.hpp collection.cpp
	$(CXX) $(CCOPTS) -DVEC$(1) -c collection.cpp -o collection$(1).o

libsim$(1).so: vecrand$(1).o interaction$(1).o constraints$(1).o collection$(1).o 
	$(CXX) $(CCOPTS) -DVEC$(1) -shared -o libsim$(1).so vecrand$(1).o interaction$(1).o constraints$(1).o collection$(1).o

endef

$(foreach target,$(VECOPTS),$(eval $(call VEC_TARGET_RULE,$(target))))

LJatoms: libsim3D.so LJatoms.cpp
	$(CXX) $(CCOPTS) -DVEC3D LJatoms.cpp -L. -lsim3D -Wl,-rpath=. -o LJatoms

LJatoms2D: libsim2D.so LJatoms.cpp
	$(CXX) $(CCOPTS) -DVEC2D LJatoms.cpp -L. -lsim2D -Wl,-rpath=. -o LJatoms2D

basicsim.zip: collection.cpp  collection.hpp  constraints.cpp  constraints.hpp  interaction.cpp  interaction.hpp  LJatoms.cpp  makefile  vec.hpp  vecrand.cpp  vecrand.hpp .vmdrc
	zip -r basicsim.zip vec.hpp {vecrand,interaction,collection,constraints}.{h,c}pp LJatoms.cpp makefile .vmdrc
