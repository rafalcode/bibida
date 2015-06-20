/* cdsck: a strict check for cds, which also outputs a filtered file */
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PDCHAR 'N' /* put this up at top, what the padding character will be */

#ifdef DBG
#define GBUF 8
#else
#define GBUF 128
#endif

#define boolean unsigned char
#define CONDREALLOC_(x, b, c, a, t); \
    if((x)>=((b)-1)) { \
        (b) += (c); \
        (a)=realloc((a), (b)*sizeof(t)); \
        memset((a)+(b)-(c), '\0', (c)*sizeof(t)); \
    }

/* another variation on above */
#define CONDREALLOC(x, b, c, a, t); \
    if((x)>=((b)-4)) { \
        (b) += (c); \
        (a)=realloc((a), (b)*sizeof(t)); \
        memset((a)+(b)-(c), '\0', (c)*sizeof(t)); \
    }

#define CONDREALLOC2(x, b, c, a1, a2, t); \
    if((x)==((b)-1)) { \
        (b) += (c); \
        (a1)=realloc((a1), (b)*sizeof(t)); \
        memset((a1)+(b)-(c), '\0', (c)*sizeof(t)); \
        (a2)=realloc((a2), (b)*sizeof(t)); \
        memset((a2)+(b)-(c), '\0', (c)*sizeof(t)); \
    }

typedef struct /* ou uov */
{
    unsigned *ua;
    unsigned ub;
    unsigned *oa;
} uo;
typedef struct /* for holding the maximums and minimums */
{ 
    unsigned mnl, mxl;
} mxn;

typedef struct /* i_s; sequence index and number of symbols */
{
    unsigned idx;
    size_t sylen; /* this is the precise symbol length of the sequence */
    char *id; /* the fasta sequence id name of the sequence */
    char *sq; /* the series of symbols that make up the sequence */
    unsigned idz, sqz;
    unsigned ibf, sbf;
} i_s; /* sequence index and number of symbols */

void usage(char *executable)
{
    printf("%s, a program to stitch a fragmented fasta file (usually a reference file).\n", executable);
    printf("Usage: %s <three arguments: 1) ref. fasta filename 2) integer, block size of each merge 3) number of merges>\n", executable);
    printf("Explanation:\n");
    printf("\tWritten because of GATK's problems dealing with reference files of very many contigs.AKA \"fragmented reference files\"\n");
    printf("\tStitching procedure described in http://gatkforums.broadinstitute.org/discussion/4774/snp-calling-using-pooled-rna-seq-datathe alignment itself.\n");
    printf("\tHere, there's a basic attempt is made to minimize number of padding N's, so the merge is done in blocks, a certain number of times.\n");
    printf("Example:\n");
    printf("\tYou have a 100k sequence file, and 90k sequences are pretty small (i.e. the fragmented ones - these numbers will be more difficult in reality).\n");
    printf("\tSo, you decide on blocks of 9k sequence, and you want 10 of those merged so that the new size will be 90010 or 9.001k.\n");
    printf("\ti.e. newsize = currentsize - (blocksize-1) * mergetimes. Command will be: \n");
    printf("\tfastitch <nmyfragmentedreffile.fa> 90000 10.\n");
    printf("\tBeware it's not easy to work out the best values of \"blocksize\" and \"mergetimes\".\n");
    printf("\tThe sequence quantity may drop, but you could get some very long sequences due to the N-padding. Use the \"fasnck\" prog to help you.\n");
    printf("Prequisite:\n");
    printf("\tYou need to know the number of sequences in your file, and have a basic clue as to how their sizes are laid out.\n");
    printf("\tYou can use the program \"fasnck\" to do this. \"fasnck\" is part of the bibida repository. Compile it with \"make fasnck\".\n");
}

int cmpsqibyl(const void *a, const void *b) /* compare sqi by length */
{
    // const int *ia = (const int *)a; // casting pointer types
    // const int *ib = (const int *)b;
    // return *ia  - *ib; /* integer comparison: returns negative if b > a and positive if a > b */
    i_s *ia = (i_s *)a;
    i_s *ib = (i_s *)b;
    // return ia->sylen  - ib->sylen; /* integer comparison: returns negative if b > a and positive if a > b: i.e. lowest values first */
    return ib->sylen  - ia->sylen; /* integer comparison: returns positive if b > a and nagetive if a > b: i.e. highest values first */
    // return (int)(100.f*ia->price - 100.f*ib->price);
    //     /* float comparison: returns negative if b > a
    //     and positive if a > b. We multiplied result by 100.0
    //     to preserve decimal fraction */
}

void prtrep(int numsq, mxn xn)
{
    printf("#SQ=%i, MXL=%u , MNL=%u\n", numsq, xn.mxl, xn.mnl);
}

void prtfaf(char *sid, char *ssq, FILE *fp) /* prints out one sequence in fasta style to a fileptr */
{
    fprintf(fp, ">");
    fprintf(fp, "%s\n", sid);
    fprintf(fp, "%s\n", ssq);
}

void gmergefirstn(i_s **sqi_, int *nsq, int n, int offset) /* gradual merge the first n squences, preparation for the progressive merge */
{
    char paddingchar=PDCHAR;
    /* want to start at the end so we can reduced number of sequences.*/
    int i;
    int numsq=*nsq;
    i_s *sqi=*sqi_;
    int noff=n+offset;
    unsigned mx=sqi[numsq-noff].sylen;
    unsigned currsz=mx*(n-1); /* total number of N's to be padded for this group */
    for(i=numsq-noff;i<numsq-offset;++i) { /* add each of the sylens to to the padded N's to work out new size of merge */
        currsz += sqi[i].sylen;
    }
#ifdef DBG
    printf("currsz= %u\n", currsz); 
#endif
    char *sqinw=malloc((1+currsz)*sizeof(char));
    memcpy(sqinw, sqi[numsq-1-offset].sq, sqi[numsq-1-offset].sylen*sizeof(char));
    char *tpos=sqinw+sqi[numsq-1-offset].sylen;
    for(i=numsq-2-offset; i>=numsq-noff;--i) {
        memset(tpos, paddingchar, mx*sizeof(char));
        tpos += mx;
        memcpy(tpos, sqi[i].sq, sqi[i].sylen*sizeof(char));
        tpos += sqi[i].sylen;
    }
    sqinw[currsz]='\0'; /* this string not null terminated yet */
#ifdef DBG
    printf("sqinwsq= %s\n", sqinw); 
#endif

    /* prepare to overwrite */
    sqi[numsq-noff].sylen=currsz;
    sqi[numsq-noff].sqz=currsz+1;
    sqi[numsq-noff].sq=realloc(sqi[numsq-noff].sq, (1+currsz)*sizeof(char));
    memcpy(sqi[numsq-noff].sq, sqinw, sqi[numsq-noff].sqz*sizeof(char));

    /* That's the first part over, now we need to move the offsetted ends of the array up accoding to the "n" entries merged
     * i.e. push all the skipped end sequences up n, so we can delete the end */
    /* TODO: this is an expensive operation, you don't have to do it each time you stitch, but, it's an optmization
     * which I want to avoid, right now */
    for(i=1; i<=offset; ++i) {
        sqi[numsq-noff+i].sylen = sqi[numsq-offset+i-1].sylen;
        sqi[numsq-noff+i].sqz = sqi[numsq-offset+i-1].sqz;
        sqi[numsq-noff+i].sq=realloc(sqi[numsq-noff+i].sq, sqi[numsq-offset+i-1].sqz*sizeof(char));
        memcpy(sqi[numsq-noff+i].sq, sqi[numsq-offset+i-1].sq, sqi[numsq-offset+i-1].sqz*sizeof(char));
    }

    /* the whole array is now going to change size: perform! */
    numsq = numsq -n +1; /* reflect new numsq value */
    for(i=numsq;i<numsq+n-1;++i) {
        free(sqi[i].id);
        free(sqi[i].sq);
    }
    sqi=realloc(sqi, numsq*sizeof(i_s));
    *sqi_=sqi; /* and, because we've realloc'd. sqi could easily have move in memory, and therefore it must be reassigned to incoming argument */
    *nsq=numsq; /* ditto numsq */

    free(sqinw);
}

void prtasf(i_s *sqi, int sz, FILE *fp) /* print all sequences to file */
{
    int i;
    for(i=0;i<sz;++i) {
        fprintf(fp, ">%s\n", sqi[i].id);
        fprintf(fp, "%s\n", sqi[i].sq);
    }
}

void uniquelens(i_s *sqi, unsigned numsq)
{
    unsigned char new;
    unsigned i, j;
    unsigned ai=0;
    uo uov;
    uov.ub=GBUF;
    uov.ua=calloc(uov.ub, sizeof(unsigned));
    uov.oa=calloc(uov.ub, sizeof(unsigned));
    for(i=0; i<numsq;++i) {
        new=1;
        for(j=0; j<=ai;++j) {
            if(uov.ua[j] == sqi[i].sylen) {
                uov.oa[j]++;
                new=0;
                break;
            }
        }
        if(new) {
            CONDREALLOC2(ai, uov.ub, GBUF, uov.ua, uov.oa, unsigned);
            uov.ua[ai]=sqi[i].sylen;
            uov.oa[ai]++;
            ai++;
        }
    }
#ifdef DBG
    printf("number of different sequence lengths: %u\n", ai);
    printf("vals: "); 
    for(j=0; j<ai;++j)
        printf("%5u ", uov.ua[j]);
    printf("\n"); 
    printf("ocs: "); 
    for(j=0; j<ai;++j)
        printf("%5u ", uov.oa[j]);
    printf("\n"); 
#endif

    free(uov.ua);
    free(uov.oa);
}

int main(int argc, char *argv[])
{
    /* argument accounting: remember argc, the number of arguments, _includes_ the executable */
    if(argc!=4) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    /* now, these  are the core of the merge op, but actually they can only be certain values */
    int blsz=atoi(argv[2]); /* this is the number of a sequences we're going to merge in one go: so call it a block */
    int mergetimes=atoi(argv[3]); /* this number of times a block merge iof size blsz will be done */

    /* general declarations */
    FILE *fin;
    char IGLINE, begline;
    size_t lidx;
    int i, c, sqidx;
    int gbuf;
    int numsq;
    i_s *sqi;
    int ididx=0;
    size_t totbases=0;
    mxn xn;
    xn.mxl=0;
    xn.mnl=0xFFFFFFFF;

    if(!(fin=fopen(argv[1], "r")) ) { /* should one check the extension of the fasta file? */
        printf("Error. Cannot open \"%s\" file.\n", argv[1]);
        exit(EXIT_FAILURE);
    }

    IGLINE=0, begline=1;
    lidx=0;

    sqidx=-1; /* this is slightly dangerous, you need very much to know what you're doing */
    gbuf=GBUF;
    sqi=malloc(gbuf*sizeof(i_s));
    for(i=gbuf-GBUF;i<gbuf;++i) {
        sqi[i].ibf=GBUF;
        sqi[i].sbf=GBUF;
        sqi[i].id=calloc(sqi[i].ibf, sizeof(char));
        sqi[i].sq=calloc(sqi[i].sbf, sizeof(char));
    }
    ididx=0;

    /* get output file ready: this will hold the stitched output file */
    size_t fnsz=1+strlen(argv[1]);
    char *tp, *tp2, *foutname=calloc(256, sizeof(char));
    char insertstr[10]="_stitched";
    foutname=realloc(foutname, (fnsz+10)*sizeof(char));
    tp=strrchr(argv[1], '.'); /* the final . position */
    tp2=strrchr(argv[1], '/'); /* the final / position, useful if files are in another directory */
    if(tp2)
        sprintf(foutname, "%.*s%s%s", (int)(tp-tp2-1), tp2+1, insertstr, tp);
    else 
        sprintf(foutname, "%.*s%s%s", (int)(tp-argv[1]), argv[1], insertstr, tp);
    FILE *fout=fopen(foutname, "w");

    while( ( (c = fgetc(fin)) != EOF) ) {
        if(c =='\n') {
            IGLINE=0;
            begline=1;
            lidx++;
        } else if( (begline==1) & (c == '>') ) { /* this condition catches the beginning of a new sequence, and uses it to prepare the nextsequence.*/
            IGLINE =1;
            begline=0; 
            if(sqidx>=0) { /* we're not interested in the first run, when sqidx=-1. Our work on the sqi array is "retrospective" */

                CONDREALLOC_(ididx, sqi[sqidx].ibf, GBUF, sqi[sqidx].id, char);
                sqi[sqidx].id[ididx]='\0';
                CONDREALLOC(sqi[sqidx].sylen, sqi[sqidx].sbf, GBUF, sqi[sqidx].sq, char);
                sqi[sqidx].sq[sqi[sqidx].sylen]='\0';
                sqi[sqidx].idz=1+ididx;
                sqi[sqidx].sqz=1+sqi[sqidx].sylen;
                totbases += sqi[sqidx].sylen;
                if(sqi[sqidx].sylen > xn.mxl)
                    xn.mxl = sqi[sqidx].sylen;
                if(sqi[sqidx].sylen < xn.mnl)
                    xn.mnl = sqi[sqidx].sylen;

                /* OK, now we have a whole sequence in the sqi[sqidx].sq variable */
                // prtfaf(sqi[sqidx].id, sqi[sqidx].sq, fout);
            }

            sqidx++;
            if(sqidx==gbuf) {
                gbuf+=GBUF;
                sqi=realloc(sqi, gbuf*sizeof(i_s));
                for(i=gbuf-GBUF;i<gbuf;++i) {
                    sqi[i].ibf=GBUF;
                    sqi[i].sbf=GBUF;
                    sqi[i].id=calloc(sqi[i].ibf, sizeof(char));
                    sqi[i].sq=calloc(sqi[i].sbf, sizeof(char));
                }
            }
            sqi[sqidx].idx=sqidx;

            /* resetting stuff */
            sqi[sqidx].sylen=0;
            ididx=0;
        } else if (IGLINE==1) {
            CONDREALLOC_(ididx, sqi[sqidx].ibf, GBUF, sqi[sqidx].id, char);
            sqi[sqidx].id[ididx++]=c;
        } else if (IGLINE==0) {
            CONDREALLOC(sqi[sqidx].sylen, sqi[sqidx].sbf, GBUF, sqi[sqidx].sq, char);
            sqi[sqidx].sq[sqi[sqidx].sylen++]=c;

        }
    }
    fclose(fin);

    /* the last sequence requires special treatment */
    CONDREALLOC_(ididx, sqi[sqidx].ibf, GBUF, sqi[sqidx].id, char);
    sqi[sqidx].id[ididx]='\0';
    CONDREALLOC_(sqi[sqidx].sylen, sqi[sqidx].sbf, GBUF, sqi[sqidx].sq, char);
    sqi[sqidx].sq[sqi[sqidx].sylen]='\0';
    sqi[sqidx].idz=1+ididx;
    sqi[sqidx].sqz=1+sqi[sqidx].sylen;
    totbases += sqi[sqidx].sylen;
    if(sqi[sqidx].sylen > xn.mxl)
        xn.mxl = sqi[sqidx].sylen;
    if(sqi[sqidx].sylen < xn.mnl)
        xn.mnl = sqi[sqidx].sylen;
    // prtfaf(sqi[sqidx].id, sqi[sqidx].sq, fout);

    /* normalize */
    numsq=sqidx+1;
    for(i=numsq;i<gbuf;++i) {
        free(sqi[i].id);
        free(sqi[i].sq);
    }
    sqi=realloc(sqi, numsq*sizeof(i_s));

    /* some checks: we've got to make sure that the blsz and the mergetimes are compatible with this numsq */
    if( (numsq-blsz*mergetimes) <0 ) {
        printf("Sorry the product of the block size (2rd arg) and mergetimes (3rd arg)\n");
        printf("is over the total number of sequences. Program can't be run. Please change these two values\n");
        printf("so that their product is either equal to, or less than, the total number of sequences.\n"); 
        exit(EXIT_FAILURE);
    }
    int nwsz=numsq-(mergetimes*(blsz-1));
    float nwzpc=100.*(float)nwsz/numsq;
    printf("INFO: The number of sequences in the new stitched file will be: %i, i.e. %3.2f%% of original.\n", nwsz, nwzpc); 
    printf("INFO: The stitched filename is now being written to your current directory, and is called \"%s\"\n", foutname);
    printf("INFO: For pretty big files, it might have taken a minute to get here. The writing out may take 10 times as long.\n");

    /* our object now is to merge the smaller sequences, so a critical first step is to sort based
     * on sequence size. Remember the struct array will be shortened, so the largest sequences should come first
     * so by simply using realloc, we can reduce the array size from its end. */
    qsort(sqi, numsq, sizeof(i_s), cmpsqibyl);

    /* OK, we going to gradually merge the sequences together */

    for(i=0;i<mergetimes;++i) {
        gmergefirstn(&sqi, &numsq, blsz, i);
#ifdef DBG
        printf("inloop idx: %i numsq val: %i\n", i, numsq); 
#endif
    }

    prtasf(sqi, numsq, fout); /* prinout to a file named stitched */
    fclose(fout);

    /* the following was because I was going to do some heuristics by calculating histograms
     * and stuff, but that's too polished ... leave until later */
    /* print report to output */
    // prtrep(numsq, xn);
    // uniquelens(sqi, numsq);

    for(i=0; i<numsq;++i) {
        free(sqi[i].id);
        free(sqi[i].sq);
    }
    free(sqi);

    free(foutname);
    return 0;
}
