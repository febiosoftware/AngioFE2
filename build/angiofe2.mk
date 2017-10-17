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

FEBIOLIBS = -Wl,--start-group $(FEBIOLIB) $(FEBIOMECH) $(FECORE)
FEBIOLIBS += $(FEBIOPLOT) $(FEBIOFLUID) $(FEBIOMIX) $(FEBIOXML) -Wl,--end-group

FEBIOLIBSO = $(FEBIOLIB) $(FEBIOMECH) $(FECORE)
FEBIOLIBSO += $(FEBIOPLOT) $(FEBIOFLUID) $(FEBIOMIX) $(FEBIOXML)

$(LIB): $(OBJ)
ifeq ($(findstring lnx,$(PLAT)),lnx)
		$(CC) $(LNKFLG) $(DBGFLG) -shared -Wl,-soname,$(SO) -o $(LIB) $(OBJ) $(FEBIOLIBS)
else ifeq ($(findstring gcc,$(PLAT)),gcc)
		$(CC) $(LNKFLG) $(DBGFLG) -shared -Wl,-soname,$(SO) -o $(LIB) $(OBJ) $(FEBIOLIBS)
else ifeq ($(findstring sky,$(PLAT)),sky)
		$(CC) $(LNKFLG) $(DBGFLG) -shared -Wl,-soname,$(SO) -o $(LIB) $(OBJ) $(FEBIOLIBS)
else
		$(CC) -dynamiclib $(FLG) -o $(LIB) $(OBJ) $(FEBIOLIBSO) $(LIBS)
endif

%.o: $(ANGDIR)AngioFE2/%.cpp
	$(CC) $(INC) $(DBGFLG) $(FLG) -MMD -c -o $@ $<

clean:
	$(RM) *.o *.d $(LIB)

-include $(DEP)
