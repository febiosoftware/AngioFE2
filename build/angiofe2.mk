SRC = $(wildcard $(ANGDIR)AngioFE2/*.cpp)
OBJ = $(patsubst $(ANGDIR)AngioFE2/%.cpp, %.o, $(SRC))
DEP = $(patsubst $(ANGDIR)AngioFE2/%.cpp, %.d, $(SRC))


SO = libangiofe2_$(PLAT).$(SFX)
LIB = $(ANGDIR)build/lib/$(SO)

FECORE = $(FEBLIB)/libfecore_$(PLAT).a

FEBIOMECH = $(FEBLIB)/libfebiomech_$(PLAT).a

FEBIOMIX = $(FEBLIB)/libfebiomix_$(PLAT).a

FEBIOLIBS = $(FEBIOMECH) $(FECORE) $(FEBIOMIX)

$(LIB): $(OBJ)
ifeq ($(findstring lnx,$(PLAT)),lnx)
		$(CC) $(LNKFLG) $(DBGFLG) -shared -Wl,-soname,$(SO) -o $(LIB) $(OBJ) $(FEBIOLIBS)
else ifeq ($(findstring gcc,$(PLAT)),gcc)
		$(CC) $(LNKFLG) $(DBGFLG) -shared -Wl,-soname,$(SO) -o $(LIB) $(OBJ) $(FEBIOLIBS)
else
		$(CC) -dynamiclib $(FLG) -o $(LIB) $(OBJ) $(FEBIOLIBS)
endif

%.o: $(ANGDIR)AngioFE2/%.cpp
	$(CC) $(INC) $(DBGFLG) $(FLG) -MMD -c -o $@ $<

clean:
	$(RM) *.o *.d $(LIB)

-include $(DEP)
