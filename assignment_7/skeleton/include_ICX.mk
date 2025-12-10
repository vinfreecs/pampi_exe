CC   = mpiicx
GCC  = gcc
LINKER = $(CC)

ifeq ($(ENABLE_OPENMP),true)
OPENMP   = -qopenmp
endif

ifeq ($(ENABLE_MPI),true)
MPI  = -D_MPI
endif

VERSION  = --version
CFLAGS   =  -O3 -xHost -qopt-zmm-usage=high -std=c99 $(OPENMP) $(MPI) -Wno-unused-command-line-argument
LFLAGS   = $(OPENMP) 
DEFINES  = -D_GNU_SOURCE -D_MPI
INCLUDES =
LIBS     =
