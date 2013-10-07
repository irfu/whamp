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
CC   = /bin/cc -Ae
AR   = ar r
BINDIR = .
# BINDIR = ~/bin

FFLAGS = -c -Wall

all: whamp whampengine whampenginetest

%.o: %.f90
	$(FTN) $(FFLAGS) $<

ALLSRC =energy.f90 rint.f90\
	rtay.f90  xsi.f90  output.f90  ryla.f90\
	rasy.f90  typin.f90 difu.f90   whamp.f90\
	whampe_engine.f90 comoutput.f90\
	input.f90 ctof.f90 av.f90 comin.f90 comcout.f90

OBJS =  energy.o  rint.o\
	rtay.o  xsi.o  output.o  ryla.o\
	rasy.o  typin.o difu.o   whamp.o\
	 input.o  ctof.o av.o comin.o comcout.o

OBJSTEST =  energy.o  rint.o\
	rtay.o  xsi.o  output.o  ryla.o\
	rasy.o  typin.o difu.o   whamp_engine_test.o\
	whamp_engine.o comoutput.o\
	 input.o  ctof.o av.o comin.o comcout.o

comcout.o: comcout.f90
	${FTN}  ${FFLAGS}  comcout.f90

comin.o: comin.f90
	${FTN}  ${FFLAGS}  comin.f90

comoutput.o: comoutput.f90
	${FTN}  ${FFLAGS}  comoutput.f90

whamp:  ${OBJS} 
	${LD}  ${OBJS} -o ${BINDIR}/whamp

whampenginetest:  ${OBJSTEST} 
	${LD}  ${OBJSTEST} -o ${BINDIR}/whamp_engine_test

whampengine: ${OBJSTEST}
	${AR} libwhamp.a whamp_engine.o comin.o comoutput.o comcout.o difu.o xsi.o energy.o ryla.o rtay.o rasy.o rint.o 

energy.o  output.o  difu.o: comcout.o
xsi.o   output.o  difu.o: comin.o
whamp_engine.o: comoutput.o
whamp_engine_test.o: comoutput.o

clean :
	rm -f *.mod *.pcl *.pc *.o $(main) *.inc *.vo *.d *~

