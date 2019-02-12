

###HOST=linux-gfortran
#HOST=macosx-gfortran
#HOST=macosx-gfortran-openmp


ifeq ($(HOST),macosx-gfortran)
  PROJECT = muller_macosx
  OBJSUF = o
  MODSUF = mod
  FC = gfortran-8
  FFLAGS = -O2
  LDFLAGS = 
  #FLINK=gfortran -o $(PROJECT) -static
  FLINK = gfortran-8 -o $(PROJECT)
endif

ifeq ($(HOST),macosx-gfortran-openmp)
  PROJECT = muller_macosx_openmp
  OBJSUF = o
  MODSUF = mod
  FC = gfortran-8
  FFLAGS = -O2 -fopenmp -w
  LDFLAGS = -fopenmp -Wl,-stack_size,0x80000000
  FLINK = gfortran-8 -o $(PROJECT)
  export OMP_NUM_THREADS=4
  export OMP_STACKSIZE=2048M
endif

ifeq ($(HOST),linux-gfortran-openmp)
  PROJECT = muller_linux_openmp
  OBJSUF = o
  MODSUF = mod
  FC = gfortran
  FFLAGS = -O2 -fopenmp -w
  LDFLAGS =
  FLINK = gfortran -fopenmp -o $(PROJECT)
  export OMP_NUM_THREADS=32
  export OMP_STACKSIZE=4096M
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




#
# SOURCE FILE LIST
#

fsrcs = test-muller.f \
  atrirouts.f \
  atritools3.f \
  patchmatc4.f \
  dotcross3d.f \
  inter3dn.f \
  rsolid.f \
  emdyadic.f \
  emplanew.f \
  emrouts2.f \
  emabrot2.f \
  dfft.f \
  hjfuns3d.f \
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
