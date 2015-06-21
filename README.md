bibida
======

The name of this repository derives from the words "BIg BIoinformatics DAta", and its objective is to both record useful notes and tips as well to script and program some procedures which deal especially with large data bionformatics directories and/or files.

## Initial orientation

The combined effect of higher speeds and lower prices in the sequencing of DNA has led to a many-fold increase in both the size and quantity of the files that are now habitually generated in bioformatics studies. Some of the most basic command-line tools that - certainly, in the case of Unix/linux - one uses when carrying out routine file handling run into trouble when handling such large and numerous files. Example of such situations are:

* Wanting to "take a quick look" at a dataset, and attempting to open it with a text editor. "Editing" is simply not what you want, or, need to do. You only need to view, to navigate. So, other tools are required. That said, the main problem is that most editors will try to read the whole file into memory.
* entering a results directory and typing "ls" only to find that the system appears to hang because the directory has over a million files.
* a variable on the above is wanting to execute a command and using tabc-ompletion to identify the file ... tab-completions appears to hang.

The upshot is that the everyday routine of file administration doesn't work with big bio data, and we need to change our approach.

# Tools in this repository

## ltar
Sometimes, in fact, your files are small, but you find yourself with a very high number of them. Surreptitiously they consume your hardisk not by eating up storage space on it, but rather by eating up the available address points, or inodes, of the file system. So, what you would really like to do, is bunch them all up into one big file. One venerably old tool for doing this is tar. But in big data cases it may take a long time. Now, there is a library that was later developed, so there's a suggestion that it could be faster, because you "tar up" at a "lower level". Could that be true? ltar is an attempt to check that. The earliest version just tars up files in the current directory, and does not enter any subdirectories.

## yafasmzr
A summarizer (yet another) for fasta files, another one of the many, though hopefully via one pass simple character scanning in C, it should be quite fast. The main benefit of this is being able to handle multifasta files, and summarizing basic information in the most compact possible way to screen.

## fasnck
Fasta file sanity checker. Basically gives a summary of the sequence sizes and the quantity of sequences in your fasta file

## faszck
Fasta file SIZE sanity check, a variation on the above.

## fastitch
If you've read my post in the GATK forums (http://gatkforums.broadinstitute.org/discussion/5709/a-reference-causes-excessive-runtime-on-genotypegvcfs#latest),
you're probably looking for this program. It needs to be compiled with GNU's c compiler on a linux system, but it has no dependencies, and might even work with Pelle's C.
Run without arguments to see the help: it basically gets fragmented fasta files and merges the smallest so reduce sequence quantity
which is what GATK is fussy, at least when it concerns the reference file.
compile with "make fastitch"
Also handy is to be able to characterize your fasta file with the fasnck, fast file sanity checker
Caveats:
* Reorders the reference file from largest to smallest sequences
* the merged sequences are then appended to the end, the first merged being the last in the new "stitched" sequence.
* simply performs the most basic of merging on the smallest contigs without caring whether they belong to each other.
* The new name of the sequence is merely one of the names of the merged contigs (extra coding can change this).
