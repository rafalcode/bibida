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
