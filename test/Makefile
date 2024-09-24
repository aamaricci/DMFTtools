TSTDIR=src
TSTBIN=bin

all: bin
	make -C $(TSTDIR)

bin:
	mkdir $(TSTBIN)
	mkdir $(TSTBIN)/DMFT_GF
# mkdir $(TSTBIN)/DMFT_WEISS_FIELD
	echo DMFT_GF > $(TSTBIN)/list_dir
# echo DMFT_WEISS_FIELD >> $(TSTBIN)/list_dir

test:
	./test.sh

clean:
	make clean -C $(TSTDIR)
	rm -rf $(TSTBIN)
