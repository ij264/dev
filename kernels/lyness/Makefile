#------------------------------------------------------------------------------
#                  Central Makefile for the program TERRA
#------------------------------------------------------------------------------

ARCH = MPIFORT #OPENGCC or MPIFORT or GFORTRAN or IFORT

# Set directory information MDLDIR    = CODE/MDLS
TARGETDIR = .
ARCHDIR   = CODE/ARCH
SRCDIR    = CODE/SRC
OBJDIR    = CODE/OBJ
MODDIR    = CODE/MDLS
LIBDIR    = CODE/LIBS
POW_SRCDIR = CODE/SRC_POWELL

# Set name executable
BIN = $(TARGETDIR)/MAKE_KERNEL

include $(ARCHDIR)/make.$(ARCH)

# Generate file set for dependency info
# Generate objects for powell inversion method tools $(OBJDIR)/commpi.o
OBJ = $(OBJDIR)/constants.o $(OBJDIR)/nrtype.o $(OBJDIR)/nr.o $(OBJDIR)/nrutil.o $(OBJDIR)/commpi.o $(OBJDIR)/sph_tools.o $(OBJDIR)/util.o $(OBJDIR)/inversion.o $(OBJDIR)/mnbrak.o $(OBJDIR)/brent.o $(OBJDIR)/linmin.o $(OBJDIR)/powell.o

OBJ += $(OBJDIR)/linpack_q.o $(OBJDIR)/ref_flds.o $(OBJDIR)/driver.o $(OBJDIR)/calc.o $(OBJDIR)/matrix.o  


# Default target is generation executable
$(BIN) : $(OBJDIR) $(MODDIR) $(LIBDIR) shtools $(OBJ)
	$(LK) -o $(BIN) $(OBJ) $(F90FLAGS) $(LFLAGS)

$(OBJDIR)/%.o : $(POW_SRCDIR)/%.f90 
	$(LK) $(F90FLAGS) -o $@ -c $< 

$(OBJDIR)/%.o : $(SRCDIR)/%.f90 
	$(LK) $(F90FLAGS) -o $@ -c $< 

# Generic rule for built process
$(OBJDIR)/%.o : $(SRCDIR)/%.f 
	$(FC) $(FFLAGS) -o $@ -c $<

shtools: ;
	cd CODE/SRC_SHTOOLS/ && make compile_mode=$(ARCH)

# Targets for housekeeping
.PHONY: clean realclean doc clean-doc

$(OBJDIR) : ;
	mkdir $(OBJDIR)

$(MODDIR) : ;
	mkdir $(MODDIR)

$(LIBDIR) : ;
	mkdir $(LIBDIR)

clean:
	$(RM) $(BIN)
	$(RM) $(MODDIR)/*.mod
	$(RM) $(LIBDIR)/*.a
	$(RM) $(OBJDIR)/*.o

realclean: clean
	$(RM) $(BIN)
	$(RM) -r $(OBJDIR)
	$(RM) -r $(MODDIR)
	$(RM) -r $(LIBDIR)
	$(RM) CODE/SRC_SHTOOLS/*.o

