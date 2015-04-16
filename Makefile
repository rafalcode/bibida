# BIBIDA Makefile
# small utilities  dealing with big DNA/PROTEIN datasets.
#
CC=gcc
DBGCFLAGS=-g -Wall# -pg # note the gprof option
CFLAGS=-O3
LIBS=-ltar

EXES=ltar ltar_dbg yafasumzr mulfaint fasnck

# ltar, code to use libtar .. in the very vain hope that it will be fast than tar itself!
ltar: ltar.c
	${CC} ${CFLAGS} -o $@ $^ $(LIBS)

# debug version.
ltar_dbg: ltar.c
	${CC} ${DBGCFLAGS} -o $@ $^ $(LIBS)

# Yet Another Fasta Summarizer, pretty much similar to emboss' infoseq
yafasumzr: yafasumzr.c
	${CC} ${CFLAGS} -o $@ $^

# Multiple fasta summarizer, same as above but print out is even more compact
# NOT WORKING ... usual probs, you refactor, you smash!
mulfaint: mulfaint.c
	${CC} ${CFLAGS} -o $@ $^

# Multple-fasta, multiple sequence fasta length histogrammer, used to be a sanity check, thereofre the "sn".
fasnck: fasnck.c
	${CC} ${CFLAGS} -o $@ $^

# Codon checker
cdsck: cdsck.c
	${CC} ${CFLAGS} -o $@ $^

.PHONY: clean

clean:
	rm -f ${EXES} vgcore.* core*
