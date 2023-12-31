##############################################################################
#
#   make all
#       Compile program in the current directory.
#
#   make clean
#       Remove the lib and module files in the LOCAL directory.
#
###############################################################################

compile_mode = gfortran

SHELL=/bin/sh

LIBTOOL = libtool
LIBTOOLFLAGS = -static
AR = ar
ARFLAGS = -r
RLIB = ranlib
RLIBFLAGS =
LDFLAGS = -s

LIBNAME = SHTOOLS
PROG = ../LIBS/lib$(LIBNAME).a

ifeq ($(compile_mode),MPIFORT)
# Default gfortran flags
F95=gfortran
F95FLAGS ?= -m64 -fPIC -O3 -ffast-math
else ifeq ($(compile_mode),OPENGCC)
# Default gfortran flags
F95=gfortran
F95FLAGS ?= -m64 -fPIC -O3 -ffast-math
else ifeq ($(compile_mode),GFORTRAN)
# Default gfortran flags
F95=gfortran
F95FLAGS ?= -m64 -fPIC -O3 -ffast-math
else ifeq ($(compile_mode),IFORT)
# Default intel fortran flags
F95=ifort
F95FLAGS ?= -m64 -fpp -free -O3 -Tf
endif
## Default Absoft Pro Fortran flags
#F95FLAGS ?= -m64 -O3 -YEXT_NAMES=LCS -YEXT_SFX=_ -fpic -speed_math=11 # -march=host
#else ifeq ($(F95),g95)
## Default g95 flags.
#F95FLAGS ?= -O3 -fno-second-underscore
#else ifeq ($(origin F95FLAGS), undefined)
#F95FLAGS = -m64 -O3

SRCS0 = SHTOOLS.f95 FFTW3.f95 CilmPlus.f95 CilmMinus.f95 ComputeDG82.f95 \
	ComputeDm.f95 DHaj.f95 djpi2.f95 BAtoHilm.f95 MakeGrid2D.f95 \
	GLQGridCoord.f95 MakeGridPoint.f95 MakeGridPointC.f95 \
	PlanetsConstants.f95 PlBar.f95 PlBar_d1.f95 PLegendre.f95 \
	PLegendre_d1.f95 PLegendreA.f95 PLegendreA_d1.f95 PlmBar.f95 \
	PlmBar_d1.f95 PlmIndex.f95 PlmSchmidt.f95 PlmSchmidt_d1.f95 \
	PlSchmidt.f95 PlSchmidt_d1.f95 PreGLQ.f95 Random.f95 \
	SHAdmitCorr.f95 SHBias.f95 SHBiasK.f95 SHConvertCoef.f95 \
	SHFindLWin.f95 SHGLQ.f95 SHLocalizedAdmitCorr.f95 SHMultiply.f95 \
	SHMultiTaperCSE.f95 SHMultiTaperSE.f95 SHPowerSpectra.f95 SHRead.f95 \
	SHReadJPL.f95 SHReturnTapersM.f95 SHRotateCoef.f95 \
	SHRotateRealCoef.f95 SphericalCapCoef.f95 \
	Wigner3j.f95 DownContFilter.f95 SHRead2.f95 MakeGeoidGrid.f95 \
	MakeCircleCoord.f95 SHMTCouplingMatrix.f95 \
	SHReturnTapers.f95 SHSjkPG.f95 PlON.f95 PlON_d1.f95 PlmON.f95 \
	PlmON_d1.f95 SHConfidence.f95 SHMagPowerSpectra.f95 SHPowerSpectraC.f95 \
	SHBiasAdmitCorr.f95 SHCimToVector.f95 YilmIndexVector.f95 \
	ComputeDMap.f95 SHReturnTapersMap.f95 Curve2Mask.f95 MakeEllipseCoord.f95 \
	CilmPlusRhoH.f95 CilmMinusRhoH.f95 BAtoHilmRhoH.f95 NormalGravity.f95 \
	SHMultiTaperMaskSE.f95 SHMultiTaperMaskCSE.f95 SHBiasKMask.f95

OBJS0 = SHTOOLS.o FFTW3.o CilmPlus.o CilmMinus.o ComputeDG82.o ComputeDm.o \
	DHaj.o djpi2.o BAtoHilm.o MakeGrid2D.o GLQGridCoord.o MakeGridPoint.o \
	MakeGridPointC.o PlanetsConstants.o PlBar.o PlBar_d1.o PLegendre.o \
	PLegendre_d1.o PLegendreA.o PLegendreA_d1.o PlmBar.o PlmBar_d1.o \
	PlmIndex.o PlmSchmidt.o PlmSchmidt_d1.o PlSchmidt.o PlSchmidt_d1.o \
	PreGLQ.o Random.o SHAdmitCorr.o SHBias.o SHBiasK.o SHConvertCoef.o \
	SHFindLWin.o SHGLQ.o SHLocalizedAdmitCorr.o \
	SHMultiply.o SHMultiTaperCSE.o SHMultiTaperSE.o SHPowerSpectra.o \
	SHRead.o SHReadJPL.o SHReturnTapersM.o SHRotateCoef.o \
	SHRotateRealCoef.o SphericalCapCoef.o Wigner3j.o \
	DownContFilter.o SHRead2.o MakeGeoidGrid.o \
	MakeCircleCoord.o SHMTCouplingMatrix.o SHReturnTapers.o SHSjkPG.o \
	PlON.o PlON_d1.o PlmON.o PlmON_d1.o \
	SHConfidence.o SHMagPowerSpectra.o \
	SHPowerSpectraC.o SHBiasAdmitCorr.o SHCilmToVector.o \
	YilmIndexVector.o ComputeDMap.o SHReturnTapersMap.o Curve2Mask.o \
	MakeEllipseCoord.o CilmPlusRhoH.o CilmMinusRhoH.o BAtoHilmRhoH.o \
	NormalGravity.o SHMultiTaperMaskSE.o SHMultiTaperMaskCSE.o SHBiasKMask.o
	
SRCSLAPACK = EigValSym.F95 EigValVecSym.F95 EigValVecSymTri.F95 \
	SHExpandLSQ.F95 SHMTDebias.F95 SHMTVarOpt.F95

OBJSLAPACK = EigValSym.o EigValVecSym.o EigValVecSymTri.o SHExpandLSQ.o \
	SHMTDebias.o SHMTVarOpt.o

SRCSFFTW = MakeGridDH.F95 MakeGridDHC.F95 MakeGridGLQ.F95 MakeGridGLQC.F95 \
	SHExpandDH.F95 SHExpandDHC.F95 SHExpandGLQ.F95 SHExpandGLQC.F95 \
	MakeGravGradGridDH.F95 MakeGravGridDH.F95 MakeMagGridDH.F95
	
OBJSFFTW = MakeGridDH.o MakeGridDHC.o MakeGridGLQ.o MakeGridGLQC.o \
	SHExpandDH.o SHExpandDHC.o SHExpandGLQ.o SHExpandGLQC.o \
	MakeGravGradGridDH.o MakeGravGridDH.o MakeMagGridDH.o

SRCS = $(SRCS0) $(SRCSLAPACK) $(SRCSFFTW)
OBJS = $(OBJS0) $(OBJSLAPACK) $(OBJSFFTW)


all: $(PROG)

$(PROG): $(OBJS)
	@echo
	@echo "--> Compilation of source files successful"
	@echo
	@rm -f $(PROG)
#	@$(LIBTOOL) $(LIBTOOLFLAGS) -o $(PROG) $(OBJS)
#	If you prefer to use libtool, uncomment the above line, and comment the two lines below (AR and RLIB)
	$(AR) $(ARFLAGS) $(PROG) $(OBJS)
	$(RLIB) $(RLIBFLAGS) $(PROG)
	@echo
	@echo "--> Creation of static library successful"
#	@rm -f $(OBJS)
#	@mv -f *.mod ../modules/
	@cp -f *.mod ../MDLS/
	@echo "--> Library and module files moved to ../LIB and ../MDLS"
	@echo "--> Archive created successfully"

.PHONY: clean clean-obs-mod clean-prog-mod all

clean: clean-obs-mod clean-prog-mod

clean-obs-mod:
	@-rm -f $(OBJS0)
	@-rm -f $(OBJSLAPACK)
	@-rm -f $(OBJSFFTW)
	@-rm -f *.mod

clean-prog-mod:
	@-rm -f $(PROG)
	@-rm -f ../MDLS/*.mod
	
%.o: %.f95
	$(F95) -c $(F95FLAGS) $<
%.o: %.F95
	$(F95) -c $(F95FLAGS) $(LAPACK_FLAGS) $(FFTW3_FLAGS) $<


CilmPlus.o: SHTOOLS.o
CilmMinus.o: SHTOOLS.o
ComputeDm.o: SHTOOLS.o
BAtoHilm.o: SHTOOLS.o
MakeGrid2D.o: SHTOOLS.o
GLQGridCoord.o: SHTOOLS.o
MakeGridGLQ.o: FFTW3.o SHTOOLS.o
MakeGridGLQ2.o: FFTW3.o SHTOOLS.o
MakeGridGLQC.o: FFTW3.o SHTOOLS.o
MakeGridGLQC2.o: FFTW3.o SHTOOLS.o
MakeGridDH.o: FFTW3.o SHTOOLS.o
MakeGridDH2.o: FFTW3.o SHTOOLS.o
MakeGridDHC.o: FFTW3.o SHTOOLS.o
MakeGridDHC2.o: FFTW3.o SHTOOLS.o
MakeGridPoint.o: SHTOOLS.o
MakeGridPointC.o: SHTOOLS.o
SHGLQ.o: SHTOOLS.o
SHAdmitCorr.o: SHTOOLS.o
SHBias.o: SHTOOLS.o
SHBiasK.o: SHTOOLS.o
SHBiasKMask.o: SHTOOLS.o
SHExpandDH.o: FFTW3.o SHTOOLS.o
SHExpandDH2.o: FFTW3.o SHTOOLS.o
SHExpandDHC.o: FFTW3.o SHTOOLS.o
SHExpandDHC2.o: FFTW3.o SHTOOLS.o
SHExpandGLQ.o: FFTW3.o SHTOOLS.o
SHExpandGLQ2.o: FFTW3.o SHTOOLS.o
SHExpandGLQC.o: FFTW3.o SHTOOLS.o
SHExpandGLQC2.o: FFTW3.o SHTOOLS.o
SHExpandLSQ.o: SHTOOLS.o
SHExpandLSQ2.o: SHTOOLS.o
SHFindLWin.o: SHTOOLS.o
SHLocalizedAdmitCorr.o: SHTOOLS.o
SHMultiply.o: SHTOOLS.o
SHMultiTaperCSE.o: SHTOOLS.o
SHMultiTaperSE.o: SHTOOLS.o
SHReturnTapersM.o: SHTOOLS.o
SHReturnTapers.o : SHTOOLS.o
SHRotateRealCoef.o: SHTOOLS.o
SHSjkPG.o: SHTOOLS.o
SphericalCapCoef.o: SHTOOLS.o
MakeGeoidGrid.o: SHTOOLS.o
SHMTVarOpt.o : SHTOOLS.o
PLegendreA.o  : SHTOOLS.o
PLegendreA_d1.o : SHTOOLS.o
PlmBar.o : SHTOOLS.o
PlmBar_d1.o : SHTOOLS.o
PlmSchmidt.o : SHTOOLS.o
PlmSchmidt_d1.o : SHTOOLS.o
PlmON.o : SHTOOLS.o
PlmON_d1.o : SHTOOLS.o
SHMTCouplingMatrix.o : SHTOOLS.o
SHMTDebias.o : SHTOOLS.o
SHMTDebias2.o : SHTOOLS.o
MakeGravGridDH.o : FFTW3.o SHTOOLS.o
MakeGravGridDH2.o : SHTOOLS.o
MakeMagGridDH.o: FFTW3.o SHTOOLS.o
SHBiasAdmitCorr.o : SHTOOLS.o
ComputeDMap.o : SHTOOLS.o
SHReturnTapersMap.o : SHTOOLS.o
CilmPlusRhoH.o : SHTOOLS.o
CilmMinusRhoH.o : SHTOOLS.o
BAtoHilmRhoH.o : SHTOOLS.o
MakeGravGradGridDH.o : FFTW3.o
MakeGravGradGridDH2.o : FFTW3.o
SHMultiTaperMaskCSE.o: SHTOOLS.o
SHMultiTaperMaskSE.o: SHTOOLS.o
