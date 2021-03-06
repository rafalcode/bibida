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

EXES=ltar ltar_dbg yafasumzr mulfaint fasnck cdsck_d cdsck faszck faszck_d fastitch fastitch_d fastitch_dd fastitch0 fastitch0_d faspli faspli_d chop1fa chop1fa_d uchop1fa uchop1fa_d faaln0 faaln0_d faaln0_dd pwamfa pwamfa_d pwamfa_dd pairupsnpa faf2snp fa_ck

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

# More general fasta checker based on faszck
fa_ck: fa_ck.c
	${CC} ${CFLAGS} -o $@ $^
fa_ck_d: fa_ck.c
	${CC} ${DBGCFLAGS} -o $@ $^

faszck_d: faszck.c
	${CC} ${DBGCFLAGS} -o $@ $^

# more verbose debug
faszck_dd: faszck.c
	${CC} ${DBG2CFLAGS} -o $@ $^

# Stats on fasta files that are actually alignments. seqret will convert for you
# this version can only take one file, which is probably going to be the most common way to use it.
faaln0: faaln0.c
	${CC} ${CFLAGS} -o $@ $^
faaln0_d: faaln0.c
	${CC} ${DBGCFLAGS} -o $@ $^
faaln0_dd: faaln0.c
	${CC} ${DBG2CFLAGS} -o $@ $^

# Stats on fasta files that are actually alignments. seqret will convert for you
# Used to be called faaln2 .. but that's a dreadful name
# pwamfa: Pairwise (verb) A MultiFasta Alignment. Note how the A is for alignment, not FA in fasta. That's important.
pwamfa: pwamfa.c
	${CC} ${CFLAGS} -o $@ $^
pwamfa_d: pwamfa.c
	${CC} ${DBGCFLAGS} -o $@ $^
pwamfa_dd: pwamfa.c
	${CC} ${DBG2CFLAGS} -o $@ $^
# Taking up  from pwamfa ... more user oriented though. 
pairupsnpa: pairupsnpa.c
	${CC} ${DBGCFLAGS} -o $@ $^

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


# Cue a supercilious creation ... almost the same as chop1fa ... but the opposite. Excludes the desired pieces.
uchop1fa: uchop1fa.c
	${CC} ${CFLAGS} -o $@ $^

uchop1fa_d: uchop1fa.c
	${CC} ${DBGCFLAGS} -o $@ $^

faf2snp: faf2snp.c
	${CC} ${CFLAGS} -o $@ $^

.PHONY: clean

clean:
	rm -f ${EXES} vgcore.* core*
