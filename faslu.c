/* cdsck: a strict check for cds, which also outputs a filtered file */
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
    printf("Usage: %s <one argument: name of fasta file>\n", executable);
    printf("Explanation:\n");
    printf("\tWritten because of GATK's problems dealing with refrence files of very many contigs.\n");
    printf("\tStitching procedure described in http://gatkforums.broadinstitute.org/discussion/4774/snp-calling-using-pooled-rna-seq-datathe alignment itself.\n");
    printf("\tHere, an attempt is made to minimize number of padding N's.\n");
}

int cmpsqibyl(const void *a, const void *b) /* compare sqi by length */
{
    // const int *ia = (const int *)a; // casting pointer types
    // const int *ib = (const int *)b;
    // return *ia  - *ib; /* integer comparison: returns negative if b > a and positive if a > b */
    i_s *ia = (i_s *)a;
    i_s *ib = (i_s *)b;
    return ia->sylen  - ib->sylen; /* integer comparison: returns negative if b > a and positive if a > b: i.e. lowest values first */
    // return (int)(100.f*ia->price - 100.f*ib->price);
    //     /* float comparison: returns negative if b > a
    //     and positive if a > b. We multiplied result by 100.0
    //     to preserve decimal fraction */
}

void prtrep(unsigned numsq, mxn xn)
{
    printf("#SQ=%u, MXL=%u , MNL=%u\n", numsq, xn.mxl, xn.mnl);
}

void prtfaf(char *sid, char *ssq, FILE *fp) /* prints out one sequence in fasta style to a fileptr */
{
    fprintf(fp, ">");
    fprintf(fp, "%s\n", sid);
    fprintf(fp, "%s\n", ssq);
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
    if(argc!=2) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    /* general declarations */
    FILE *fin;
    char IGLINE, begline;
    size_t lidx;
    int i, c, sqidx;
    int gbuf;
    unsigned numsq;
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
    qsort(sqi, numsq, sizeof(i_s), cmpsqibyl);
    prtasf(sqi, numsq, fout);

    fclose(fout);

    /* print report to output */
    // prtrep(numsq, xn);
    uniquelens(sqi, numsq);

    for(i=0; i<numsq;++i) {
        free(sqi[i].id);
        free(sqi[i].sq);
    }
    free(sqi);

    free(foutname);
    return 0;
}
