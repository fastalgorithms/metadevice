

###HOST=linux-gfortran
#HOST=osx-gfortran
#HOST=osx-gfortran-openmp

NAME = fmps

ifeq ($(HOST),osx-gfortran)
  PROJECT = $(NAME)_osx
  OBJSUF = o
  MODSUF = mod
  FC = gfortran-8
  FFLAGS = -O2
  LDFLAGS = 
  #FLINK=gfortran -o $(PROJECT) -static
  FLINK = gfortran-8 -o $(PROJECT)
endif

ifeq ($(HOST),osx-gfortran-openmp)
  PROJECT = $(NAME)_osx_openmp
  OBJSUF = o
  MODSUF = mod
  FC = gfortran-8
  FFLAGS = -O2 -fopenmp -w
  LDFLAGS = -fopenmp -Wl,-stack_size,0x80000000 -framework accelerate
  FLINK = gfortran-8 -o $(PROJECT)
  export OMP_NUM_THREADS=1
  export OMP_STACKSIZE=2048M
endif

ifeq ($(HOST),osx-intel)
  PROJECT = $(NAME)_osx
  FC = ifort -c -w
  FFLAGS = -O2
  LDFLAGS = 
  FLINK = ifort -w -mkl -o $(PROJECT)
endif

ifeq ($(HOST),osx-intel-openmp)
  PROJECT = $(NAME)_osx_openmp
  FC = ifort -c -w -qopenmp
  FFLAGS = -O2
  FLINK = ifort -w -mkl=parallel -qopenmp \
    -Wl,-stack_size,0x80000000 -o $(PROJECT)
  export OMP_NUM_THREADS=4
  export OMP_STACKSIZE=2048M
endif

ifeq ($(HOST),linux-gcc-openmp)
  PROJECT = $(NAME)_linux_openmp
  OBJSUF = o
  MODSUF = mod
  FC = gfortran
  FFLAGS = -O2 -fopenmp -w
  LDFLAGS =
  FLINK = gfortran -fopenmp -o $(PROJECT) -lblas -llapack
  #ulimit -s unlimited
  #export OMP_NUM_THREADS=32
  #export OMP_STACKSIZE=4096M
endif

ifeq ($(HOST),linux-gcc)
  PROJECT = $(NAME)_linux
  OBJSUF = o
  FC = gfortran
  FFLAGS = -O2 -w
  LDFLAGS =
  FLINK = gfortran -o $(PROJECT) -lblas -llapack
endif

ifeq ($(HOST),linux-intel-openmp)
  PROJECT = $(NAME)_linux_openmp
  OBJSUF = o
  FC = ifort
  FFLAGS = -O2 -qopenmp -w
  LDFLAGS =
  FLINK = ifort -qopenmp -o $(PROJECT) -mkl
  export OMP_NUM_THREADS=32
  export OMP_STACKSIZE=4096M
endif

ifeq ($(HOST),linux-intel)
  PROJECT = $(NAME)_linux
  OBJSUF = o
  FC = ifort
  FFLAGS = -O2 -w
  LDFLAGS =
  FLINK = ifort -o $(PROJECT) -mkl
endif




#
# SOURCE FILE LIST
#
SRCDIR = ../src

f90srcs = fmps_sphere.f90 \
  $(SRCDIR)/xtri_plot.f90 \
  $(SRCDIR)/lapack_wrap.f90 \
  $(SRCDIR)/patchmatc.f90

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
  $(SRCDIR)/selfquad_new.f \
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
