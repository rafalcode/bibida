# On a linux system, type "make spot" to compile

CC=gcc
CFLAGS=-Wall
EXES=facg2wig batchfa2wig bafa2wig filtyex
INSTPATH=/usr/local/bin

# Just takes a single  fasta and coverts it to wig.
facg2wig: c5ar1.c
	${CC} ${CFLAGS} -o $@ $^
	sudo cp $@ $(INSTPATH)

# Especially aimed at very many fasta files all in one directory
batchfa2wig: c5ar2.c
	${CC} ${CFLAGS} -o $@ $^
	sudo cp $@ $(INSTPATH)


# this actually renders all the fasta files as GC data into one big wig file called All.wig and then also a all.sizes file with name and sequence size
bafa2wig: c5ar3.c
	${CC} ${CFLAGS} -o $@ $^
	sudo cp $@ $(INSTPATH)

filtyex: filtyex.c
	${CC} ${CFLAGS} -o $@ $^
	sudo cp $@ $(INSTPATH)

.PHONY: clean

clean:
	rm -f ${EXES}
