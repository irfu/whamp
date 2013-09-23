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
BINDIR = .
# BINDIR = ~/bin

FFLAGS = -c -Wall

all: whamp 

%.o: %.f90
	$(FTN) $(FFLAGS) $<

ALLSRC =energy.f90 rint.f90  wrfi.f90  blockdata.f90\
	rtay.f90  xsi.f90  output.f90  ryla.f90\
	rasy.f90  typin.f90 difu.f90   whamp.f90\
	input.f90 ctof.f90 av.f90 comin.f90 comcout.f90

OBJS =  energy.o  rint.o  wrfi.o  blockdata.o\
	rtay.o  xsi.o  output.o  ryla.o\
	rasy.o  typin.o difu.o   whamp.o\
	 input.o  ctof.o av.o comin.o comcout.o

comcout.o: comcout.f90
	${FTN}  ${FFLAGS}  comcout.f90

comin.o: comin.f90
	${FTN}  ${FFLAGS}  comin.f90

whamp:  ${OBJS} 
	${LD}  ${OBJS} -o ${BINDIR}/whamp

energy.o  wrfi.o  output.o  difu.o: comcout.o
wrfi.o    xsi.o   output.o  difu.o: comin.o

clean :
	rm -f *.mod *.pcl *.pc *.o $(main) *.inc *.vo *.d *~

