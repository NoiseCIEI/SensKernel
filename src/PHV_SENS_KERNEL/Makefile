#
INST_DIR = ../../bin
BIN  = PHV_SENS_KERNEL 

fflags =  -Wall -ffixed-line-length-none

FFLAGS = $(DBG) $(fflags)

FC = gfortran

DBG =

FOBJS =   PHV_SENS_KERNEL.o  

$(BIN) : $(FOBJS)
	$(FC)  $(FFLAGS) $(FOBJS) -o $(BIN) $(LDLIBS)

install :: $(BIN)
	install -s $(BIN) $(INST_DIR)

clean ::
	rm -f $(BIN) core $(FOBJS)
