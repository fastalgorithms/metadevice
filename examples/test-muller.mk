

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
  LDFLAGS = -fopenmp
  FLINK = gfortran -o $(PROJECT)
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
SRCDIR = ../src

f90srcs = test-muller.f90

fsrcs = $(SRCDIR)/emutils.f \
  $(SRCDIR)/atrirouts.f \
  $(SRCDIR)/atritools3.f \
  $(SRCDIR)/patchmatc4.f \
  $(SRCDIR)/dotcross3d.f \
  $(SRCDIR)/inter3dn.f \
  $(SRCDIR)/rsolid.f \
  $(SRCDIR)/emdyadic.f \
  $(SRCDIR)/emplanew.f \
  $(SRCDIR)/emrouts2.f \
  $(SRCDIR)/emabrot2.f \
  $(SRCDIR)/dfft.f \
  $(SRCDIR)/hjfuns3d.f \
  $(SRCDIR)/rotviarecur3.f \
  $(SRCDIR)/yrecursion.f \
  $(SRCDIR)/xrecursion.f \
  $(SRCDIR)/triaadap.f \
  $(SRCDIR)/ctriaadap.f \
  $(SRCDIR)/triagauc.f \
  $(SRCDIR)/triasymq.f \
  $(SRCDIR)/selfquad.f \
  $(SRCDIR)/radial.f \
  $(SRCDIR)/print.f \
  $(SRCDIR)/legendre.f \
  $(SRCDIR)/c8triadam.f \
  $(SRCDIR)/c9triadam.f \
  $(SRCDIR)/c28triadam.f \
  $(SRCDIR)/c29triadam.f \
  $(SRCDIR)/cgmres_rel.f \
  $(SRCDIR)/cgmressq_rel.f \
  $(SRCDIR)/cbicgstab_rel.f \
  $(SRCDIR)/cqrsolve.f \
  $(SRCDIR)/legeexps.f \
  $(SRCDIR)/ortho2eva.f \
  $(SRCDIR)/ortho2exps4.f \
  $(SRCDIR)/orthom.f \
  $(SRCDIR)/prini.f \
  $(SRCDIR)/prinm.f \
  $(SRCDIR)/xprini.f

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
