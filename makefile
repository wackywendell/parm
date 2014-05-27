UNAME := $(shell uname)
#CXX=${CXX}
SWIG=swig -Wextra -shadow -python -py3 -c++
CCOPTS=-Wall -O2 -fPIC #-std=c++11

INC=`python3-config --includes`

.PHONY: all 2d 3d 2dlong 3dlong printout clean wraps

all: 2d 2dlong 3d 3dlong
	@echo "making all."

2d: lib/_sim2d.so

3d: lib/_sim3d.so

2dlong: lib/_sim2dlong.so

3dlong: lib/_sim3dlong.so

wraps: lib/sim_wrap2d.cxx lib/sim_wrap2dlong.cxx lib/sim_wrap3d.cxx lib/sim_wrap3dlong.cxx

printout:
	@echo Running on \"$(UNAME)\"

clean:
	rm -f bin/LJatoms bin/LJatom2d bin/hardspheres
	cd lib; rm -f *.o *.so *.gch sim_wrap*.cxx

lib/sim_wrap2d.cxx: src/sim.i src/collection.hpp src/constraints.hpp src/interaction.hpp src/trackers.hpp src/box.hpp src/vecrand.hpp src/collection.cpp src/constraints.cpp src/interaction.cpp src/trackers.cpp src/box.cpp src/vecrand.cpp src/vec.hpp
	cd src ; $(SWIG) -DVEC2D sim.i
	mv src/sim_wrap.cxx lib/sim_wrap2d.cxx

lib/sim_wrap3d.cxx: src/sim.i src/collection.hpp src/constraints.hpp src/interaction.hpp src/trackers.hpp src/box.hpp src/vecrand.hpp src/collection.cpp src/constraints.cpp src/interaction.cpp src/trackers.cpp src/box.cpp src/vecrand.cpp src/vec.hpp
	cd src ; $(SWIG) -DVEC3D sim.i
	mv src/sim_wrap.cxx lib/sim_wrap3d.cxx

lib/sim_wrap2d.o: lib/sim_wrap2d.cxx
	$(CXX) $(CCOPTS) -DVEC2D -I src/ -c lib/sim_wrap2d.cxx -o lib/sim_wrap2d.o $(INC)

lib/sim_wrap3d.o: lib/sim_wrap3d.cxx
	$(CXX) $(CCOPTS) -DVEC3D -I src/ -c lib/sim_wrap3d.cxx -o lib/sim_wrap3d.o $(INC)

lib/_sim2d.so: lib/sim_wrap2d.o
	$(CXX) $(CCOPTS) -DVEC2D -shared lib/sim_wrap2d.o -o lib/_sim2d.so $(LIB)
	
lib/_sim3d.so: lib/sim_wrap3d.o
	$(CXX) $(CCOPTS) -DVEC3D -shared lib/sim_wrap3d.o -o lib/_sim3d.so $(LIB)

lib/sim_wrap2dlong.cxx: src/sim.i src/collection.hpp src/constraints.hpp src/interaction.hpp src/trackers.hpp src/box.hpp src/vecrand.hpp src/collection.cpp src/constraints.cpp src/interaction.cpp src/trackers.cpp src/box.cpp src/vecrand.cpp src/vec.hpp
	cd src ; $(SWIG) -DVEC2D -DLONGFLOAT sim.i
	mv src/sim_wrap.cxx lib/sim_wrap2dlong.cxx

lib/sim_wrap3dlong.cxx: src/sim.i src/collection.hpp src/constraints.hpp src/interaction.hpp src/trackers.hpp src/box.hpp src/vecrand.hpp src/collection.cpp src/constraints.cpp src/interaction.cpp src/trackers.cpp src/box.cpp src/vecrand.cpp src/vec.hpp
	cd src ; $(SWIG) -DVEC3D -DLONGFLOAT sim.i
	mv src/sim_wrap.cxx lib/sim_wrap3dlong.cxx

lib/sim_wrap2dlong.o: lib/sim_wrap2dlong.cxx
	$(CXX) $(CCOPTS)  -Wconversion -DVEC2D -DLONGFLOAT -I src/ -c lib/sim_wrap2dlong.cxx  -o lib/sim_wrap2dlong.o $(INC)

lib/sim_wrap3dlong.o: lib/sim_wrap3dlong.cxx
	$(CXX) $(CCOPTS) -Wconversion -DVEC3D -DLONGFLOAT -I src/ -c lib/sim_wrap3dlong.cxx  -o lib/sim_wrap3dlong.o $(INC)

lib/_sim2dlong.so: lib/sim_wrap2dlong.o
	$(CXX) $(CCOPTS) -Wconversion -DVEC2D -DLONGFLOAT -shared lib/sim_wrap2dlong.o -o lib/_sim2dlong.so $(LIB)
	
lib/_sim3dlong.so: lib/sim_wrap3dlong.o
	$(CXX) $(CCOPTS) -Wconversion -DVEC3D -DLONGFLOAT -shared lib/sim_wrap3dlong.o -o lib/_sim3dlong.so $(LIB)

VECOPTS := 2D 3D

define VEC_TARGET_RULE
lib/vecrand$(1).o: src/vec.hpp src/vecrand.hpp src/vecrand.cpp
	$(CXX) $(CCOPTS) -DVEC$(1) -c src/vecrand.cpp -o lib/vecrand$(1).o

lib/box$(1).o: src/vec.hpp src/vecrand.hpp src/box.hpp src/box.cpp
	$(CXX) $(CCOPTS) -DVEC$(1) -c src/box.cpp -o lib/box$(1).o

lib/trackers$(1).o: src/vec.hpp src/vecrand.hpp src/box.hpp src/trackers.hpp src/trackers.cpp
	$(CXX) $(CCOPTS) -DVEC$(1) -c src/trackers.cpp -o lib/trackers$(1).o

lib/interaction$(1).o: src/vec.hpp src/vecrand.hpp src/box.hpp src/trackers.hpp src/interaction.hpp src/interaction.cpp
	$(CXX) $(CCOPTS) -DVEC$(1) -c src/interaction.cpp -o lib/interaction$(1).o

lib/constraints$(1).o: src/vec.hpp src/vecrand.hpp src/box.hpp src/trackers.hpp src/interaction.hpp src/constraints.hpp src/constraints.cpp
	$(CXX) $(CCOPTS) -DVEC$(1) -c src/constraints.cpp -o lib/constraints$(1).o

lib/collection$(1).o: src/vec.hpp src/vecrand.hpp src/box.hpp src/trackers.hpp src/interaction.hpp src/collection.hpp src/constraints.hpp src/collection.cpp
	$(CXX) $(CCOPTS) -DVEC$(1) -c src/collection.cpp -o lib/collection$(1).o

lib/libsim$(1).so: lib/vecrand$(1).o lib/box$(1).o lib/trackers$(1).o lib/interaction$(1).o lib/constraints$(1).o lib/collection$(1).o 
	$(CXX) $(CCOPTS) -DVEC$(1) -shared -o lib/libsim$(1).so lib/box$(1).o lib/trackers$(1).o lib/vecrand$(1).o lib/interaction$(1).o lib/constraints$(1).o lib/collection$(1).o

endef

$(foreach target,$(VECOPTS),$(eval $(call VEC_TARGET_RULE,$(target))))

bin/LJatoms: lib/libsim3D.so src/LJatoms.cpp
	$(CXX) $(CCOPTS) -DVEC3D src/LJatoms.cpp -Llib -lsim3D -Wl,-rpath=. -o bin/LJatoms

bin/LJatoms2d: lib/libsim2D.so src/LJatoms.cpp
	$(CXX) $(CCOPTS) -DVEC2D src/LJatoms.cpp -Llib -lsim2D -Wl,-rpath=. -o bin/LJatoms2d

bin/hardspheres: lib/libsim3D.so src/hardspheres.cpp
	$(CXX) $(CCOPTS) -DVEC3D src/hardspheres.cpp -Llib -lsim3D -Wl,-rpath=. -o bin/hardspheres
