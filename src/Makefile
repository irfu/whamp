#*
# Makefile     -- Makefile to build the WHAMP program.
#
# Created 910916 by Kjell Ronnmark 
#
#===============================================================================
#
#
FTN  = gfortran
LD   = gfortran
AR   = ar r
BINDIR = .

FFLAGS = -c -Wall -fPIC

all: whamp libwhamp.a whamp_engine_test

# General rule (most f90 built separatly)
%.o: %.f90
	$(FTN) $(FFLAGS) $<


# Special dependencies
energy.o  output.o  difu.o: comcout.o
xsi.o     output.o  difu.o: comin.o
whamp_engine.o: comoutput.o
whamp_engine_test.o: comoutput.o


OBJS =  energy.o rint.o\
	rtay.o xsi.o output.o ryla.o\
	rasy.o typin.o difu.o\
	whamp.o\
	input.o ctof.o av.o comin.o comcout.o

OBJSTEST = energy.o rint.o\
	rtay.o xsi.o output.o ryla.o\
	rasy.o typin.o difu.o\
	whamp_engine_test.o\
	whamp_engine.o comoutput.o\
	input.o ctof.o av.o comin.o comcout.o

OBJSLIB = whamp_engine.o comin.o comoutput.o\
	comcout.o difu.o xsi.o energy.o\
	ryla.o rtay.o rasy.o rint.o


whamp:  ${OBJS}
	${LD} ${OBJS} -o ${BINDIR}/$@

whamp_engine_test: ${OBJSTEST}
	${LD} ${OBJSTEST} -o ${BINDIR}/$@

libwhamp.a: ${OBJSLIB}
	${AR} $@ ${OBJSLIB}


# Clean up everything
clean :
	rm -f *.mod *.pcl *.pc *.o *.inc *.vo *.d *~ libwhamp.a whamp whamp_engine_test

