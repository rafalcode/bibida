/* This is multifasta summarizer aimed at being fast, and also giving out only very compressed summary stat
 *
 * - It tolerates empty lines and will not count them
 * - It only reports DNA, and tolerates big caps and small caps symbols.
 * - Ambiguity symbols are dealt with ..., for exampel Y increments the thw coutn for th two pyramidines, etc.
 * - Any other symbols will view as anomalous and will call your sequence "AnoSQ" instead of "SQ", but these symbols will be excluded from the total count.
 * - It reports the sequence indexin the file, its total base count, and then the foru percetnages of the different nucleotides.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GBUF 4
#define SBUF 8
#define SSZ 5 /* DNA symbol alphabet size is 4 ? Of course, but the fifth is for N's. Ambiguous, anomalous, none */

typedef struct /* the i_s struct */
{
    char *n;
    unsigned nsz;
    unsigned int idx;
    size_t tsz; /* because summing over the ambiguity characters will include the double and tripled nature of some of them */
    size_t sy[SSZ];
} i_s; /* sequence index and number of symbols */

void prti_s(i_s *sqisz, int sz)
{
    int i;
    char *sqgood;

    size_t tsz;
    for(i=0;i<sz;++i) {
        if(sqisz[i].sy[SSZ-1] != 0)
            sqgood="AnoSQ";
        else
            sqgood="SQ";
        tsz=sqisz[i].tsz;

            printf(">%s (sz:%u) ", sqisz[i].n, sqisz[i].nsz);
            printf("| %s#%i=TOT: %zu A:%.3f C:%.3f G:%.3f T:%.3f ?:%zu\n", sqgood, i, tsz, (float)sqisz[i].sy[0]/tsz, (float)sqisz[i].sy[1]/tsz, (float)sqisz[i].sy[2]/tsz, (float)sqisz[i].sy[3]/tsz, sqisz[i].sy[4]);
    }
    printf("\n"); 
}

int main(int argc, char *argv[])
{
    /* argument accounting: remember argc, the number of arguments, _includes_ the executable */
    if(argc!=2) {
        printf("Error. Pls supply 1 argument:; fasta file name.\n");
        exit(EXIT_FAILURE);
    }
    FILE *fin;
    if(!(fin=fopen(argv[1], "r")) ) {
        printf("Error. Cannot open the presented filename\n");
        exit(EXIT_FAILURE);
    }

    char IGLINE=0, begline=1;
    int nidx=0;
    size_t lidx=0;
    int i, c, sqidx=-1; /* this is slightly dangerous, you need very much to knwo what you're doing */
    int gbuf=GBUF;
    int sbuf=SBUF;
    i_s *sqisz=malloc(gbuf*sizeof(i_s));
    for(i=0;i<gbuf;++i) 
        sqisz[i].n=calloc(sbuf, sizeof(char));

    while( ( (c = fgetc(fin)) != EOF) ) {
        if(c =='\n') {
            if(IGLINE==1) {
                sqisz[sqidx].n[nidx]='\0';
                sqisz[sqidx].n=realloc(sqisz[sqidx].n, (nidx+1)*sizeof(char));
                sbuf=SBUF;
                sqisz[sqidx].nsz=nidx;
                nidx=0;
            }
            IGLINE=0;
            begline=1;
            lidx++;
        } else if( (begline==1) & (c == '>') ) {
            IGLINE =1;
            begline=0; 
            sqidx++;
            if(sqidx==gbuf) {
                gbuf+=GBUF;
                sqisz=realloc(sqisz, gbuf*sizeof(i_s));
                for(i=gbuf-GBUF;i<gbuf;++i)
                    sqisz[i].n=calloc(sbuf, sizeof(char));
            }
            sqisz[sqidx].idx=sqidx;
            sqisz[sqidx].tsz=0;
            for(i=0;i<SSZ;++i)
                sqisz[sqidx].sy[i]=0;
        } else if ((IGLINE==1) & (begline==0)) {
            if(nidx==sbuf-1) {
                sbuf+=SBUF;
                sqisz[sqidx].n=realloc(sqisz[sqidx].n, sbuf*sizeof(char));
                memset(sqisz[sqidx].n+sbuf-SBUF, 0, SBUF*sizeof(char));
            }
            sqisz[sqidx].n[nidx]=c;
            nidx++;
        } else if (IGLINE==0) {
            sqisz[sqidx].tsz++;
            switch(c) {
                case 'A': case 'a':
                    sqisz[sqidx].sy[0]++; break;
                case 'C': case 'c':
                    sqisz[sqidx].sy[1]++; break;
                case 'G': case 'g':
                    sqisz[sqidx].sy[2]++; break;
                case 'T': case 't':
                    sqisz[sqidx].sy[3]++; break;
                case 'M': case 'm': 
                    sqisz[sqidx].sy[0]++;
                    sqisz[sqidx].sy[1]++; break;
                case 'R': case 'r': 
                    sqisz[sqidx].sy[0]++;
                    sqisz[sqidx].sy[2]++; break;
                case 'W': case 'w':
                    sqisz[sqidx].sy[0]++;
                    sqisz[sqidx].sy[3]++; break;
                case 'S': case 's':
                    sqisz[sqidx].sy[1]++;
                    sqisz[sqidx].sy[2]++; break;
                case 'Y': case 'y':
                    sqisz[sqidx].sy[1]++;
                    sqisz[sqidx].sy[3]++; break;
                case 'K': case 'k':
                    sqisz[sqidx].sy[2]++;
                    sqisz[sqidx].sy[3]++; break;
                case 'V': case 'v':
                    sqisz[sqidx].sy[0]++;
                    sqisz[sqidx].sy[1]++;
                    sqisz[sqidx].sy[2]++; break;
                case 'H': case 'h':
                    sqisz[sqidx].sy[0]++;
                    sqisz[sqidx].sy[1]++;
                    sqisz[sqidx].sy[3]++; break;
                case 'D': case 'd': 
                    sqisz[sqidx].sy[0]++;
                    sqisz[sqidx].sy[2]++;
                    sqisz[sqidx].sy[3]++; break;
                case 'B': case 'b':
                    sqisz[sqidx].sy[1]++;
                    sqisz[sqidx].sy[2]++;
                    sqisz[sqidx].sy[3]++; break;
                case 'N': case 'n':
                    sqisz[sqidx].sy[0]++;
                    sqisz[sqidx].sy[1]++;
                    sqisz[sqidx].sy[2]++;
                    sqisz[sqidx].sy[3]++; break;
                default:
                    sqisz[sqidx].tsz--; break;
            }
        }
    }
    fclose(fin);
    int numsq=sqidx+1, numano=0;
    for(i=0;i<numsq;++i) 
        if(sqisz[i].sy[SSZ-1])
            numano++;
    for(i=numsq;i<gbuf;++i) 
        free(sqisz[i].n);
    sqisz=realloc(sqisz, numsq*sizeof(i_s));
    prti_s(sqisz, numsq);
    /* the summary comes at the end because otherwise, with many sequences, it goes off-screen */
    printf("Lines in file: %zu. Number of sequences: %i, of which %i have non-ACGTN symbols.\n", lidx, numsq, numano);

    for(i=0;i<numsq;++i) 
        free(sqisz[i].n);
    free(sqisz);
    return 0;
}
