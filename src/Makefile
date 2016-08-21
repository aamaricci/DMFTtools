include ../make.inc

REV = $(shell git rev-parse HEAD)

DIR = ./dmft_tools_portable

.SUFFIXES:.f90 .f

all: version dmftt lib


dmftt: DMFT_TOOLS_VERSION.o DMFT_FFTGF.o DMFT_FFTAUX.o DMFT_MISC.o DMFT_CONVERGENCE.o DMFT_VECTORS.o DMFT_TIGHT_BINDING.o DMFT_CTRL_VARS.o DMFT_GLOC.o DMFT_EKIN.o DMFT_TOOLS.o


lib:
	@rsync -avP *.mod $(INC_TARGET)/
	ar crv $(LIB_DMFTT) `ls *.o | sort | uniq` 
	ranlib $(LIB_DMFTT)
	@echo 'DMFT_Tools library is done.........'
	@echo ' '

version:
	@echo "Actual version is: $(REV)"
	@echo 'character(len=41),parameter,public :: sf_version = "$(REV)"' > dmft_tools_version.inc

clean:
	@echo ""
	@echo "Cleaning ..."
	rm -fv *.o *.mod
	@echo ""
	@echo "done."
	@echo ""

portable:
	./create_portable.sh $(DIR)


cmdexit:
	exit 1



.f90.o:
	$(FC) -c $(FCOPT) -I$(INC_SCIFOR)/ -I./ $(MOPT)./ $< -o ./$@

.f.o:
	$(FC) -c $(FCOPT)  $< -o ./$@
