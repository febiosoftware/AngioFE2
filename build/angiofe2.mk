SRC = $(wildcard $(ANGDIR)AngioFE2/*.cpp)
OBJ = $(patsubst $(ANGDIR)AngioFE2/%.cpp, %.o, $(SRC))
DEP = $(patsubst $(ANGDIR)AngioFE2/%.cpp, %.d, $(SRC))


SO = libangiofe2_$(PLAT).$(SFX)
LIB = $(ANGDIR)build/lib/$(SO)

FECORE = $(FEBLIB)/libfecore_$(PLAT).a

FEBIOLIB = $(FEBLIB)/libfebiolib_$(PLAT).a

FEBIOMECH = $(FEBLIB)/libfebiomech_$(PLAT).a

FEBIOPLOT = $(FEBLIB)/libfebioplot_$(PLAT).a

FEBIOFLUID = $(FEBLIB)/libfebiofluid_$(PLAT).a

FEBIOMIX = $(FEBLIB)/libfebiomix_$(PLAT).a

FEBIOXML = $(FEBLIB)/libfebioxml_$(PLAT).a

INTELROOT = $(subst /mkl,,$(MKLROOT))/compiler
INTEL_INC = $(INTELROOT)/include
INTEL_LIB = $(INTELROOT)/lib/intel64

MKL_PATH = $(MKLROOT)/lib/intel64
MKL_LIB = -Wl,--start-group $(MKL_PATH)/libmkl_intel_lp64.a
MKL_LIB += $(MKL_PATH)/libmkl_intel_thread.a $(MKL_PATH)/libmkl_core.a -Wl,--end-group
MKL_LIB += -liomp5 -pthread -lz
#MKL_LIB += $(INTEL_LIB)/libiomp5.a -pthread -lm -ldl

FEBIOLIBS = -Wl,--start-group $(FEBIOLIB) $(FEBIOMECH) $(FECORE)
FEBIOLIBS += $(FEBIOPLOT) $(FEBIOFLUID) $(FEBIOMIX) $(FEBIOXML) $(MKL_LIB)

FEBIOLIBSO = $(FEBIOLIB) $(FEBIOMECH) $(FECORE)
FEBIOLIBSO += $(FEBIOPLOT) $(FEBIOFLUID) $(FEBIOMIX) $(FEBIOXML)

$(LIB): $(OBJ)
ifeq ($(findstring lnx,$(PLAT)),lnx)
		$(CC) $(FLG) $(LNKFLG) $(DBGFLG) -shared -Wl,-soname,$(SO) -o $(LIB) $(OBJ) $(FEBIOLIBS)
else ifeq ($(findstring gcc,$(PLAT)),gcc)
		$(CC) $(FLG) $(LNKFLG) $(DBGFLG) -shared -Wl,-soname,$(SO) -o $(LIB) $(OBJ) $(FEBIOLIBS)
else ifeq ($(findstring sky,$(PLAT)),sky)
		$(CC) $(FLG) $(LNKFLG) $(DBGFLG) -shared -Wl,-soname,$(SO) -o $(LIB) $(OBJ) $(FEBIOLIBS)
else
		$(CC) -dynamiclib $(FLG) -o $(LIB) $(OBJ) $(FEBIOLIBSO) $(LIBS)
endif

%.o: $(ANGDIR)AngioFE2/%.cpp
	$(CC) $(INC) $(DBGFLG) $(FLG) -MMD -c -o $@ $<

clean:
	$(RM) *.o *.d $(LIB)

-include $(DEP)
