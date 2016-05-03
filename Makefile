# BIBIDA Makefile
# small utilities  dealing with big DNA/PROTEIN datasets.
# Beware wthat many of the first few program do not store the sequence in memory, and only store statistics on the sequence.
# Use fastitch as template if you also want the sequence.
#
CC=gcc
DBGCFLAGS=-g -Wall -DDBG
DBG2CFLAGS=-g -Wall -DDBG2
CFLAGS=-O3
LIBS=-ltar

EXES=ltar ltar_dbg yafasumzr mulfaint fasnck cdsck_d cdsck faszck faszck_d fastitch fastitch_d fastitch_dd fastitch0 fastitch0_d faspli chop1fa

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

# Multple-fasta, multiple sequence fasta length size checker. THis is very similar to the histogram
# but actually outputs the unique sequence lengthsof the fasta togther with how often that
# length occurred.
faszck: faszck.c
	${CC} ${CFLAGS} -o $@ $^

faszck_d: faszck.c
	${CC} ${DBGCFLAGS} -o $@ $^

# more verbose debug
faszck_dd: faszck.c
	${CC} ${DBG2CFLAGS} -o $@ $^

# Codon checker
cdsck: cdsck.c
	${CC} ${CFLAGS} -o $@ $^

# Codon checker
cdsck_d: cdsck.c
	${CC} ${DBGCFLAGS} -o $@ $^

# fastitch: takes bits from cdsck ot create a prog
# Note for the git logs, this used to be calle faslu i.e. slurp fasta, during devel
fastitch: fastitch.c
	${CC} ${CFLAGS} -o $@ $^

fastitch_d: fastitch.c
	${CC} ${DBGCFLAGS} -o $@ $^

# Split a multifasta file into a 
# what still needs to be done .. niceer naming
# dealing with DOS terminated newlines.
faspli: faspli.c
	${CC} ${CFLAGS} -o $@ $^

faspli2_d: faspli2.c
	${CC} ${DBGCFLAGS} -o $@ $^

fastitch0: fastitch2.c
	${CC} ${CFLAGS} -o $@ $^

fastitch0_d: fastitch2.c
	${CC} ${DBGCFLAGS} -o $@ $^

chop1fa: chop1fa.c
	${CC} ${CFLAGS} -o $@ $^

chop1fa2_d: chop1fa2.c
	${CC} ${DBGCFLAGS} -o $@ $^

.PHONY: clean

clean:
	rm -f ${EXES} vgcore.* core*
