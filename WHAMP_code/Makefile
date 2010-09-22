#*
# Makefile     -- Makefile to build the WHAMP program.
#
# Created 910916 by Kjell Ronnmark 
#
#===============================================================================
#
#
FTN  = g77
LD   = g77
CC   = /bin/cc -Ae
BINDIR = .
# BINDIR = ~/bin

FFLAGS = -g -c -fno-f90 -fugly-complex
.f.o:
	${FTN}  ${FFLAGS}  $<

ALLSRC =energy.f  rint.f  wrfi.f  blockdata.f\
	rtay.f  xsi.f  output.f  ryla.f\
	rasy.f  typin.f difu.f   whamp.f\
	comcout.h  comin.h input.f ctof.f av.f

OBJS =  energy.o  rint.o  wrfi.o  blockdata.o\
	rtay.o  xsi.o  output.o  ryla.o\
	rasy.o  typin.o difu.o   whamp.o\
	 input.o  ctof.o av.o

whamp:  ${OBJS}
	${LD}  ${OBJS} -o ${BINDIR}/whamp

energy.o  wrfi.o  output.o  difu.o: comcout.h
wrfi.o    xsi.o   output.o  difu.o: comin.h

