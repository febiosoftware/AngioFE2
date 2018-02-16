ANGDIR = $(dir $(CURDIR))

# This needs to be set in order to find the febio .h files.
FEBINC = $(FEBIO_HOME)
# This needs to be set in order to find the febio library files.
FEBLIB = $(FEBIO_HOME)/build/lib

export CC = icpc

lnx64: PLAT = lnx64
lnx64d: PLAT = lnx64d
lnx64s: PLAT = lnx64s
lnx32: PLAT = lnx32
lnx32d: PLAT = lnx32d
osx: PLAT = osx
osxs: PLAT = osxs
gcc: PLAT = gcc
sky: PLAT = sky
lnx64clean: PLAT = lnx64
lnx64dclean: PLAT = lnx64d
lnx64sclean: PLAT = lnx64s
lnx32clean: PLAT = lnx32
lnx32dclean: PLAT = dlnx32
osxclean: PLAT = osx
osxsclean: PLAT = osxs
gccclean: PLAT = gcc
skyclean: PLAT = sky

gcc: export CC = g++
sky: export CC = g++	

export PLAT
export FEBLIB
export ANGDIR

# Assumes the Intel C++ compiler.  Use g++ for the gnu compiler and run "make gcc".


export FLG = -O3 -fPIC -std=c++11

# Can be removed if not using the Intel compiler.
export LNKFLG = -lz -static-intel -no-intel-extensions -lz

export INC = -I$(FEBINC)
export DBGFLG = -DNDEBUG
export LIBS = -lz

gcc: export LNKFLG = -lz
sky: export LNKFLG = -lz

lnx64 lnx32 gcc sky:
	( export DBGFLG=-DNDEBUG; cd $(PLAT)/AngioFE2; $(MAKE) -f ../../angiofe2.mk SFX="so" )

lnx64d lnx64s lnx32d:
	( cd $(PLAT)/AngioFE2; $(MAKE) -f ../../angiofe2.mk SFX="so" )

osx osxs:
	( cd $(PLAT); $(MAKE) -f ../angiofe2.mk SFX="dylib" )

lnx64clean lnx32clean gccclean lnx64dclean lnx32dclean skyclean:
	( cd $(PLAT)/AngioFE2; $(MAKE) -f ../../angiofe2.mk SFX="so" clean )

lnx64sclean:
	( cd $(PLAT); $(MAKE) -f ../angiofe2.mk SFX="so" clean )

osxclean osxsclean:
	( cd $(PLAT); $(MAKE) -f ../angiofe2.mk SFX="dylib" clean )

.PHONY: lnx64 lnx32 osx osxs gcc lnx64d lnx64s lnx32d 
