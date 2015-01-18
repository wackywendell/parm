UNAME := $(shell uname)
#CXX=${CXX}
SWIG=swig -Wextra -shadow -python -py3 -c++
CCOPTS=-I src -Wall -O2 -fPIC -Wconversion -Wno-sign-conversion #-std=c++11

INC=`python3-config --includes`

.PHONY: all 2d 3d 2dlong 3dlong printout clean wraps py

all: bin/LJatoms bin/LJatoms2d bin/hardspheres
	@echo "making all."

py: 2d 3d 2dlong 3dlong

2d: pyparm/_sim2d.so

3d: pyparm/_sim3d.so

2dlong: pyparm/_sim2dlong.so

3dlong: pyparm/_sim3dlong.so

wraps: pyparm/sim_wrap2d.cxx pyparm/sim_wrap2dlong.cxx pyparm/sim_wrap3d.cxx pyparm/sim_wrap3dlong.cxx

printout:
	@echo Running on \"$(UNAME)\"

clean:
	rm -f bin/LJatoms bin/LJatoms2d bin/hardspheres
	cd src; rm -f *.o *.so *.gch sim_wrap*.cxx
	cd lib; rm -f *.o *.so *.gch sim_wrap*.cxx
	cd pyparm; rm -f *.o *.so *.gch sim_wrap*.cxx
	rm -f src/sim.py

lib:
	mkdir -p lib

bin:
	mkdir -p bin

pyparm/sim_wrap2d.cxx: src/sim.i src/collection.hpp src/constraints.hpp src/interaction.hpp src/trackers.hpp src/box.hpp src/vecrand.hpp src/collection.cpp src/constraints.cpp src/interaction.cpp src/trackers.cpp src/box.cpp src/vecrand.cpp src/vec.hpp
	cd src ; $(SWIG) -DVEC2D sim.i
	mv src/sim_wrap.cxx pyparm/sim_wrap2d.cxx
	mv src/sim2d.py pyparm/d2.py

pyparm/sim_wrap3d.cxx: src/sim.i src/collection.hpp src/constraints.hpp src/interaction.hpp src/trackers.hpp src/box.hpp src/vecrand.hpp src/collection.cpp src/constraints.cpp src/interaction.cpp src/trackers.cpp src/box.cpp src/vecrand.cpp src/vec.hpp
	cd src ; $(SWIG) -DVEC3D sim.i
	mv src/sim_wrap.cxx pyparm/sim_wrap3d.cxx
	mv src/sim3d.py pyparm/d3.py

pyparm/sim_wrap2d.o: pyparm/sim_wrap2d.cxx
	$(CXX) $(CCOPTS) -DVEC2D -I src/ -c pyparm/sim_wrap2d.cxx -o pyparm/sim_wrap2d.o $(INC)

pyparm/sim_wrap3d.o: pyparm/sim_wrap3d.cxx
	$(CXX) $(CCOPTS) -DVEC3D -I src/ -c pyparm/sim_wrap3d.cxx -o pyparm/sim_wrap3d.o $(INC)

pyparm/_sim2d.so: pyparm/sim_wrap2d.o
	$(CXX) $(CCOPTS) -DVEC2D -shared pyparm/sim_wrap2d.o -o pyparm/_sim2d.so $(LIB)

pyparm/_sim3d.so: pyparm/sim_wrap3d.o
	$(CXX) $(CCOPTS) -DVEC3D -shared pyparm/sim_wrap3d.o -o pyparm/_sim3d.so $(LIB)

pyparm/sim_wrap2dlong.cxx: src/sim.i src/collection.hpp src/constraints.hpp src/interaction.hpp src/trackers.hpp src/box.hpp src/vecrand.hpp src/collection.cpp src/constraints.cpp src/interaction.cpp src/trackers.cpp src/box.cpp src/vecrand.cpp src/vec.hpp
	cd src ; $(SWIG) -DVEC2D -DLONGFLOAT sim.i
	mv src/sim_wrap.cxx pyparm/sim_wrap2dlong.cxx
	mv src/sim2dlong.py pyparm/d2long.py

pyparm/sim_wrap3dlong.cxx: src/sim.i src/collection.hpp src/constraints.hpp src/interaction.hpp src/trackers.hpp src/box.hpp src/vecrand.hpp src/collection.cpp src/constraints.cpp src/interaction.cpp src/trackers.cpp src/box.cpp src/vecrand.cpp src/vec.hpp
	cd src ; $(SWIG) -DVEC3D -DLONGFLOAT sim.i
	mv src/sim_wrap.cxx pyparm/sim_wrap3dlong.cxx
	mv src/sim3dlong.py pyparm/d3long.py

pyparm/sim_wrap2dlong.o: pyparm/sim_wrap2dlong.cxx
	$(CXX) $(CCOPTS)  -Wconversion -DVEC2D -DLONGFLOAT -I src/ -c pyparm/sim_wrap2dlong.cxx  -o pyparm/sim_wrap2dlong.o $(INC)

pyparm/sim_wrap3dlong.o: pyparm/sim_wrap3dlong.cxx
	$(CXX) $(CCOPTS) -Wconversion -DVEC3D -DLONGFLOAT -I src/ -c pyparm/sim_wrap3dlong.cxx  -o pyparm/sim_wrap3dlong.o $(INC)

pyparm/_sim2dlong.so: pyparm/sim_wrap2dlong.o
	$(CXX) $(CCOPTS) -Wconversion -DVEC2D -DLONGFLOAT -shared pyparm/sim_wrap2dlong.o -o pyparm/_sim2dlong.so $(LIB)

pyparm/_sim3dlong.so: pyparm/sim_wrap3dlong.o
	$(CXX) $(CCOPTS) -Wconversion -DVEC3D -DLONGFLOAT -shared pyparm/sim_wrap3dlong.o -o pyparm/_sim3dlong.so $(LIB)

VECOPTS := 2D 3D

define VEC_TARGET_RULE
lib/vecrand$(1).o: lib src/vec.hpp src/vecrand.hpp src/vecrand.cpp
	$(CXX) $(CCOPTS) -DVEC$(1) -c src/vecrand.cpp -o lib/vecrand$(1).o

lib/box$(1).o: lib src/vec.hpp src/vecrand.hpp src/box.hpp src/box.cpp
	$(CXX) $(CCOPTS) -DVEC$(1) -c src/box.cpp -o lib/box$(1).o

lib/trackers$(1).o: lib src/vec.hpp src/vecrand.hpp src/box.hpp src/trackers.hpp src/trackers.cpp
	$(CXX) $(CCOPTS) -DVEC$(1) -c src/trackers.cpp -o lib/trackers$(1).o

lib/interaction$(1).o: lib src/vec.hpp src/vecrand.hpp src/box.hpp src/trackers.hpp src/interaction.hpp src/interaction.cpp
	$(CXX) $(CCOPTS) -DVEC$(1) -c src/interaction.cpp -o lib/interaction$(1).o

lib/constraints$(1).o: lib src/vec.hpp src/vecrand.hpp src/box.hpp src/trackers.hpp src/interaction.hpp src/constraints.hpp src/constraints.cpp
	$(CXX) $(CCOPTS) -DVEC$(1) -c src/constraints.cpp -o lib/constraints$(1).o

lib/collection$(1).o: lib src/vec.hpp src/vecrand.hpp src/box.hpp src/trackers.hpp src/interaction.hpp src/collection.hpp src/constraints.hpp src/collection.cpp
	$(CXX) $(CCOPTS) -DVEC$(1) -c src/collection.cpp -o lib/collection$(1).o

lib/libsim$(1).so: lib lib/vecrand$(1).o lib/box$(1).o lib/trackers$(1).o lib/interaction$(1).o lib/constraints$(1).o lib/collection$(1).o
	$(CXX) $(CCOPTS) -DVEC$(1) -shared -o lib/libsim$(1).so lib/box$(1).o lib/trackers$(1).o lib/vecrand$(1).o lib/interaction$(1).o lib/constraints$(1).o lib/collection$(1).o

endef

$(foreach target,$(VECOPTS),$(eval $(call VEC_TARGET_RULE,$(target))))

bin/LJatoms: bin lib/libsim3D.so src/bin/LJatoms.cpp
	$(CXX) $(CCOPTS) -DVEC3D src/bin/LJatoms.cpp -Llib -lsim3D -Wl,-rpath "lib" -o bin/LJatoms

bin/LJatoms2d: bin lib/libsim2D.so src/bin/LJatoms.cpp
	$(CXX) $(CCOPTS) -DVEC2D src/bin/LJatoms.cpp -Llib -lsim2D -Wl,-rpath "lib" -o bin/LJatoms2d

bin/hardspheres: bin lib/libsim3D.so src/bin/hardspheres.cpp
	$(CXX) $(CCOPTS) -DVEC3D src/bin/hardspheres.cpp -Llib -lsim3D -Wl,-rpath "lib" -o bin/hardspheres

bin/hardspheres2: bin lib/libsim3D.so src/bin/hardspheres2.cpp
	$(CXX) $(CCOPTS) -DVEC3D src/bin/hardspheres2.cpp -Llib -lsim3D -Wl,-rpath "lib" -o bin/hardspheres2

bin/hardspheres3: bin lib/libsim3D.so src/bin/hardspheres3.cpp
	$(CXX) $(CCOPTS) -DVEC3D src/bin/hardspheres3.cpp -Llib -lsim3D -Wl,-rpath "lib" -o bin/hardspheres3
