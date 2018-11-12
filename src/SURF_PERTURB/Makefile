#
INST_DIR = ../../bin
BIN=SURF_PERTURB   

fflags =  -Wall -ffixed-line-length-none

FFLAGS = $(DBG) $(fflags)

FC = gfortran

DBG =

FOBJS =   SURF_PERTURB.o  surfa.o  calcul_deep.o init_deep.o flat1.o

$(BIN) : $(FOBJS)
	$(FC)  $(FFLAGS) $(FOBJS) -o $(BIN) $(LDLIBS)

install :: $(BIN)
	install -s $(BIN) $(INST_DIR)

clean ::
	rm -f $(BIN) core $(FOBJS)
