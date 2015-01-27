# BIBIDA Makefile
# small utilities  dealing with big DNA/PROTEIN datasets.
#
CC=gcc
DBGCFLAGS=-g -Wall# -pg # note the gprof option
CFLAGS=-O3
LIBS=-ltar

EXES=ltar ltar_dbg jafasumzr

# ltar, code to use libtar .. in the very vain hope that it will be fast than tar itself!
ltar: ltar.c
	${CC} ${CFLAGS} -o $@ $^ $(LIBS)

# debug version.
ltar_dbg: ltar.c
	${CC} ${DBGCFLAGS} -o $@ $^ $(LIBS)

# Just Another Fasta Summarizer
jafasumzr: jafasumzr.c
	${CC} ${CFLAGS} -o $@ $^


.PHONY: clean

clean:
	rm -f ${EXES} vgcore.* core*
