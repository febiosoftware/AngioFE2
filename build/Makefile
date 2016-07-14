ANGDIR = $(dir $(CURDIR))

# This needs to be set in order to find the febio .h files.
FEBINC = $(FEBIO_HOME)
# This needs to be set in order to find the febio library files.
FEBLIB = $(FEBIO_HOME)/build/lib

lnx64: PLAT = lnx64
lnx64d: PLAT = lnx64d
lnx32: PLAT = lnx32
lnx32d: PLAT = lnx32d
osx: PLAT = osx
gcc: PLAT = gcc
lnx64clean: PLAT = lnx64
lnx32clean: PLAT = lnx32
osxclean: PLAT = osx
gccclean: PLAT = gcc

export PLAT
export FEBLIB
export ANGDIR

# Assumes the Intel C++ compiler.  Use g++ for the gnu compiler and run "make gcc".
export CC = icpc
export FLG = -O3 -fPIC -stdc++11

# Can be removed if not using the Intel compiler.
export LNKFLG = -static-intel -no-intel-extensions
export INC = -I$(FEBINC) -I$(ANGDIR)AngioLib
export DBGFLG = 

lnx64 lnx32 gcc:
	( export DBGFLG=-DNDEBUG; cd $(PLAT)/AngioFE2; $(MAKE) -f ../../angiofe2.mk SFX="so" )

lnx64d lnx32d:
	( cd $(PLAT)/AngioFE2; $(MAKE) -f ../../angiofe2.mk SFX="so" )

osx:
	( cd $(PLAT)/AngioFE2; $(MAKE) -f ../../angiofe2.mk SFX="dylib" )

lnx64clean lnx32clean gccclean:
	( cd $(PLAT)/AngioFE2; $(MAKE) -f ../../angiofe2.mk SFX="so" clean )

osxclean:
	( cd $(PLAT)/AngioFE2; $(MAKE) -f ../../angiofe2.mk SFX="dylib" clean )

.PHONY: lnx64 lnx32 osx gcc
