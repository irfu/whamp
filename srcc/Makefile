#*
# Makefile     -- Makefile to build the WHAMP C translation.
#
# Mirrors src/Makefile: same three targets (whamp, libwhamp.a,
# whamp_engine_test) with the same object lists.
#
#===============================================================================

CC   = gcc
AR   = ar r
BINDIR = .

CFLAGS = -c -Wall -O2

all: whamp libwhamp.a whamp_engine_test

# General rule
%.o: %.c whamp.h
	$(CC) $(CFLAGS) $<


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
	${CC} ${OBJS} -lm -o ${BINDIR}/$@

whamp_engine_test: ${OBJSTEST}
	${CC} ${OBJSTEST} -lm -o ${BINDIR}/$@

libwhamp.a: ${OBJSLIB}
	${AR} $@ ${OBJSLIB}


# Clean up everything
clean :
	rm -f *.o libwhamp.a whamp whamp_engine_test
