FC=mpif90
FFLAG = -O2 -ffree-line-length-none
FPPMPI=-cpp -D_MPI
FFLAGS= ${FFLAG} ${FPPMPI}

TSTBIN=../../bin/DMFT_GF

GLOB_INC:=$(shell pkg-config --cflags dmft_tools scifor)
GLOB_LIB:=$(shell pkg-config --libs dmft_tools scifor | sed  "s/;/ /g"  | sed 's/\\/  /g' )


OBJS=../ASSERTING.o
DT_INIT_SOURCES=$(wildcard *.f90)

DT_INIT_EXES=$(DT_INIT_SOURCES:.f90=)

all: $(DT_INIT_EXES)

$(DT_INIT_EXES):$(DT_INIT_SOURCES)
	$(FC) $(FFLAGS) $(OBJS) $@.f90  -o $(TSTBIN)/$@.x -I../ $(GLOB_INC) $(GLOB_LIB)

clean:
	rm -rf *.mod *.o 
