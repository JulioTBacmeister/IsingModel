##################################
#  $Id: GNUmakefile,v 1.24 2008/11/30 02:42:01 bacmj Exp $
##################################

#EXEDIR=.
#EXEN=clstr.x

SRCS=           		\
	iblock.F90 \
	flipspin.F90 \
	energy_ccc.F90 \
	diags.F90 \
	renorm_utils.F90 \
	ising_4x4_partfun.F90 \
	ising_NxN_partfun.F90 \
	ising_MC_partfun.F90 \
	ising_MC_rlzn.F90 \
	rg_NxN_rlzn.F90 \
	rg_NxNp_rlzn.F90 \
	rg_4x4_rlzn.F90 \
	rg_MC_rlzn.F90 \
	rg_MCb_rlzn.F90 \
	rng_sub2.F90 \
	linsys.F90 \
	rng_drv.F90 \
	ising_NxN_drv.F90 \
	ising_MC_drv.F90 \
	rg_4x4_drv.F90 \
	rg_MC_drv.F90 \
	minrg_drv.F90


COMPILER=gfortran
HOST=$(shell hostname)
UNAME=$(shell uname -s)

OBJS= ${SRCS:.F90=.o}

F77_FLAGS=

STD_FLAGS=-c

ifeq ($(COMPILER),gfortran)
  #STD_FLAGS := $(STD_FLAGS) -O3 -fdollar-ok -DGFORTRAN
  STD_FLAGS := $(STD_FLAGS) -g -fdollar-ok -fbounds-check -fdefault-real-8 -DGFORTRAN
endif
ifeq ($(COMPILER),pgf95)
  STD_FLAGS := $(STD_FLAGS) -g
endif

F90_FLAGS=$(PFLAGS) $(STD_FLAGS) $(USER_FDEFS)

CPP_FLAGS=-P -DALPHA_MACH

ifeq ($(UNAME),Linux)
  LAFLAGS=-L/usr/lib64 -llapack -lblas
endif
ifeq ($(UNAME),Darwin)
  LAFLAGS=-L/Users/juliob/lapack-3.8.0 -llapack -lrefblas
endif


#-----------------------------------------------------
# These files (RHS) are used by almost everyone.
# So, safer just to rebuild everything if they change.
#-----------------------------------------------------
%.o: %.F90
	$(COMPILER) $(F90_FLAGS) $<

objtest: GNUmakefile
	echo $(OBJS)

rng: $(OBJS)
	$(COMPILER) iblock.o rng_sub2.o rng_drv.o -o rng.x

isinxn: $(OBJS)
	$(COMPILER) ising_NxN_partfun.o ising_nxn_drv.o -o isinxn.x

isimc: $(OBJS)
	$(COMPILER) energy_ccc.o ising_MC_partfun.o ising_MC_rlzn.o ising_mc_drv.o -o isimc.x

ising: $(OBJS)
	$(COMPILER) diags.o ising_MC_rlzn.o ising_MC_drv.o -o ising.x

rgmc: $(OBJS)
	$(COMPILER) $(LAFLAGS) diags.o rg_MCb_rlzn.o rg_NxN_rlzn.o rg_4x4_rlzn.o rg_MC_drv.o -o rgmc.x

rg4x4: $(OBJS)
	$(COMPILER) $(LAFLAGS) renorm_utils.o rg_MCb_rlzn.o rg_NxNp_rlzn.o rg_4x4_rlzn.o rg_4x4_drv.o -o rg4x4.x

minrg: $(OBJS)
	$(COMPILER) $(LAFLAGS) renorm_utils.o minrg_drv.o -o minrg.x

linsys: $(OBJS)
	$(COMPILER) $(LAFLAGS) linsys.o -o linsys.x

testlapack:
	$(COMPILER) -fdefault-real-8 testlapack.F90 -L/Users/juliob/lapack-3.8.0 -llapack -lrefblas  -o testlapack.x


cleancode : 
	/bin/rm -f *.o *.a *.x *.mod *.dvi *.ps *.pdf *.ofl *.aux *.log *~

clean : 
	/bin/rm -f *.o *.a *.x *.mod *.dvi *.ps *.pdf *.ofl *.aux *.log *~ fort.* share/*~ idlpros/*~

enviro : 
	@echo "HOST = "$(HOST)
	@echo "OS   = "$(UNAME)
	@echo "LAPACK flags = "$(LAFLAGS)
