# In short: this is very close to supporting Python 2, but supporting 2 and 3 at
# the same time was too much for me and my `makefile`.

# Most of this library is Python-version agnostic, except for the `makefile`,
# which includes `-py3` as a swig option by default. I did not want to remove
# this, as it adds useful functionality, but I also could not find a good way to
# include it only when Python 3 is wanted.

# I tried using a setup.py that invokes swig, but that was even more difficult.

UNAME := $(shell uname)
#CXX=${CXX}
SWIG=swig -Wextra -shadow -python -py3 -c++
CCOPTS=-I src -Wall -O2 -fPIC -Wconversion -Wno-sign-conversion -std=c++98
BINOPTS:=-Llib

INC=`python3-config --includes`

.PHONY: all py2d py3d py2dlong py3dlong printout clean wraps py doc ghp

all: bin/LJatoms2d bin/LJatoms3dlong bin/hardspheres
	@echo "making all."

py: py2d py3d py2dlong py3dlong

wraps: pyparm/sim_wrap2d.cxx pyparm/sim_wrap2dlong.cxx pyparm/sim_wrap3d.cxx pyparm/sim_wrap3dlong.cxx

printout:
	@echo Running on \"$(UNAME)\"
	
ghp: doc
	echo 'parm.lostinmyterminal.com' > doc/html/CNAME
	ghp-import doc/html
	
doc: doc/html/index.html
	
doc/html/index.html: src/*.cpp src/*.hpp src/*.md src/bin/*.cpp pyparm/examples/*.py Doxyfile README.md
	doxygen Doxyfile >/dev/null
	
clean:
	rm -f bin/* lib/*
	cd src; rm -f *.o *.so *.gch sim_wrap*.cxx
	cd pyparm; rm -f *.o *.so *.gch sim_wrap*.cxx
	rm -f src/sim.py

lib:
	mkdir -p lib

bin:
	mkdir -p bin

# We now define the four options: (2d or 3d) + (doubles or long doubles)
# SFX is the suffix for object files, and MODNAME is the module name for the python modules
# OPTSET is the appropriate option set for each version

VECOPTS := 2 3

FLOATOPTS := long notlong

define TARGET_RULES

$(eval NDIM=$(1))
$(if $(findstring notlong,$(2)), $(eval FLT=), $(eval FLT=long))
$(if $(findstring notlong,$(2)), $(eval FLTOPT=), $(eval FLTOPT=-DLONGFLOAT))

$(eval OPTSET=-DVEC$(NDIM)D $(FLTOPT))
$(eval SFX:=$(NDIM)d$(FLT))
$(eval MODNAME:=d$(NDIM)$(FLT))

#-------------------------------------------------------------------------------
# The python modules

py$(SFX): pyparm/_sim$(SFX).so

# We use -DSWIG_TYPE_TABLE=sim$(SFX) so that types from sim2d and sim3d with the
# same names do not end up treated as the same types. This can cause issues with
# e.g. AtomVec, if both sim2d and sim3d are imported at the same time.
# There isn't any need for the two to interoperate, so we keep them separate.
pyparm/sim_wrap$(SFX).cxx: src/swig_header.h src/sim.i src/array.i src/collection.hpp src/constraints.hpp src/interaction.hpp src/trackers.hpp src/box.hpp src/vecrand.hpp src/collection.cpp src/constraints.cpp src/interaction.cpp src/trackers.cpp src/box.cpp src/vecrand.cpp
	cd src ; $(SWIG) $(OPTSET) -DSWIG_TYPE_TABLE=sim$(SFX) sim.i
	(cat src/swig_header.h ; echo ; echo ; cat src/sim_wrap.cxx) > $$@
	rm src/sim_wrap.cxx
	mv src/sim$(SFX).py pyparm/$(MODNAME).py

pyparm/sim_wrap$(SFX).o: pyparm/sim_wrap$(SFX).cxx
	$(CXX) $(CCOPTS) $(OPTSET) -DSWIG_TYPE_TABLE=sim$(SFX) $(INC) -I src/ -c $$< -o $$@

pyparm/_sim$(SFX).so: pyparm/sim_wrap$(SFX).o
	$(CXX) $(CCOPTS) $(OPTSET) -shared $$< -o $$@ $(LIB)

#-------------------------------------------------------------------------------
# The C++ modules
lib/vecrand$(SFX).o: src/vecrand.cpp src/vecrand.hpp | lib
	$(CXX) $(CCOPTS) $(OPTSET) -c $$< -o $$@

lib/box$(SFX).o: src/box.cpp src/vecrand.hpp src/box.hpp | lib
	$(CXX) $(CCOPTS) $(OPTSET) -c $$< -o $$@

lib/trackers$(SFX).o: src/trackers.cpp src/vecrand.hpp src/box.hpp src/trackers.hpp | lib
	$(CXX) $(CCOPTS) $(OPTSET) -c $$< -o $$@

lib/interaction$(SFX).o: src/interaction.cpp src/vecrand.hpp src/box.hpp src/trackers.hpp src/interaction.hpp | lib
	$(CXX) $(CCOPTS) $(OPTSET) -c $$< -o $$@

lib/constraints$(SFX).o: src/constraints.cpp src/vecrand.hpp src/box.hpp src/trackers.hpp src/interaction.hpp src/constraints.hpp | lib
	$(CXX) $(CCOPTS) $(OPTSET) -c $$< -o $$@

lib/collection$(SFX).o: src/collection.cpp src/vecrand.hpp src/box.hpp src/trackers.hpp src/interaction.hpp src/collection.hpp src/constraints.hpp | lib
	$(CXX) $(CCOPTS) $(OPTSET) -c $$< -o $$@

lib/libsim$(SFX).so: lib/vecrand$(SFX).o lib/box$(SFX).o lib/trackers$(SFX).o lib/interaction$(SFX).o lib/constraints$(SFX).o lib/collection$(SFX).o | lib
	$(CXX) $(CCOPTS) $(OPTSET) -shared $$^ -o $$@

lib/libsim$(SFX).a: lib/vecrand$(SFX).o lib/box$(SFX).o lib/trackers$(SFX).o lib/interaction$(SFX).o lib/constraints$(SFX).o lib/collection$(SFX).o | lib
	ar -r $$@ $$^

bin/LJatoms$(SFX): src/bin/LJatoms.cpp lib/libsim$(SFX).a | bin
	$(CXX) $(CCOPTS) $(OPTSET) $(BINOPTS) $$^ -o $$@

bin/packer$(SFX): src/bin/packer.cpp lib/libsim$(SFX).a | bin
	$(CXX) $(CCOPTS) $(OPTSET) $(BINOPTS) $$^ -o $$@

endef

$(foreach target1,$(VECOPTS), $(foreach target2,$(FLOATOPTS),$(eval $(call TARGET_RULES,$(target1),$(target2)))))

bin/hardspheres: src/bin/hardspheres.cpp lib/libsim3d.a | bin
	$(CXX) $(CCOPTS) $(BINOPTS) -DVEC3D $^ -o $@

bin/hardspheres2: src/bin/hardspheres2.cpp lib/libsim3d.a | bin
	$(CXX) $(CCOPTS) $(BINOPTS) -DVEC3D $^ -o $@
