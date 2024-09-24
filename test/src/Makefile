FC=mpif90
FFLAG = -O2 -ffree-line-length-none
FPPMPI=-cpp -D_MPI
FFLAGS= ${FFLAG} ${FPPMPI}

.SUFFIXES: .f90
.PHONY: clean

GLOB_INC:=$(shell pkg-config --cflags scifor)

OBJS=ASSERTING.o

all: $(OBJS) gf

gf:
	make -C ./DMFT_GF


.f90.o:	
	$(FC) $(FFLAGS) -c $< ${GLOB_INC}

clean:
	rm -rf *.mod *.o 
