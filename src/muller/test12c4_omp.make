
###HOST=linux-gfortran
HOST=macosx-gfortran


ifeq ($(HOST),macosx-gfortran)
  PROJECT = muller_macosx
  OBJSUF = o
  MODSUF = mod
  FC = gfortran
  FFLAGS = -O2
  LDFLAGS = 
  #FLINK=gfortran -o $(PROJECT) -static
  FLINK = gfortran -o $(PROJECT)
endif

# ifeq ($(HOST),linux-gfortran-openmp)

# # buggy compiler, proceed anyway, use gfortran > 4.4.0
# PROJECT=muller_lnx64_omp
# OBJSUF=o
# MODSUF=mod
# FC=gfortran -c 
# FFLAGS=-O3 --openmp
# FLINK=gfortran -o $(PROJECT) --openmp
# ### export OMP_NUM_THREADS=4
# ### export OMP_STACKSIZE=1024M

# else



f90srcs = ../../utils/lapack_wrap.f90
fsrcs = ../../utils/hank103.f \
  $(IDDIR)/dfft.f \
  $(IDDIR)/id_rand.f \
  $(IDDIR)/id_rtrans.f \



# SOURCE FILE LIST
#
vpath %.f .:../emlib:../lib:../../H3DLibraries



FSRCS = test12c4.f atritools3.f \
        patchmatc4.f dotcross3d.f inter3dn.f rsolid.f atrirouts.f \
           emdyadic.f emplanew.f \
           emrouts2.f emabrot2.f dfft.f hjfuns3d.f \
           rotviarecur3.f yrecursion.f xrecursion.f \
           triaadap.f ctriaadap.f \
           triagauc.f triasymq.f \
           selfquad.f radial.f print.f legendre.f \
           c8triadam.f \
           c9triadam.f \
           c28triadam.f \
           c29triadam.f \
           cgmres_rel.f cgmressq_rel.f cbicgstab_rel.f \
           cqrsolve.f \
           legeexps.f \
           ortho2eva.f ortho2exps4.f orthom.f \
           prini.f prinm.f xprini.f


ifeq ($(WITH_SECOND),1) 
FSRCS += second-r8.f
endif



objs = $(fsrcs:.f=.o) $(f90srcs:.f90=.o)


# rule to generate a dep file by using the C preprocessor
# (see man cpp for details on the -MM and -MT options)


.PHONY: all
all: $(objs)
	rm -f $(PROJECT)
	$(FLINK) $^ $(LDFLAGS)
	./$(PROJECT)

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(objs)
	rm -f $(PROJECT)
