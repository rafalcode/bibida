/* mulfaint.c only takes a single sequence nucleotide fasta file,
 * [reject RNA; protein and mutlifasta files]
 * converts the symbols into integers: ambiguous taken note of */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SBUF 4
#define GBUF 4
#define SSZ 2 /* CG count, first, AT count second, third are the anomalous characters */
#define HISTBUCKETSZ 25

typedef struct /* i_s; sequence index and number of symbols */
{
    unsigned int idx;
    size_t tsz; /* because summing over the ambiguity characters will include the double and tripled nature of some of them */
    size_t sy[SSZ];
    float cgp;
    char *intsq; /* the sequence symbols converted to ints (well, 1-16) */
    unsigned intsqbuf;
    unsigned ambano[2]; /* number of ambiguous symbots (first), then number of anomalous symbols */
} i_s; /* sequence index and number of symbols */

void prthist(char *histname, int *bucketarr, int numbuckets)
{
    int i;
    printf("Hist for %s: ", histname); 
    printf("lowest<--"); 
    for(i=0;i<numbuckets;++i) 
        printf("| %i ", bucketarr[i]);
    printf("|-->highest\n"); 
}

int *hist_cg(i_s *sqisz, int sz, float mxcg, float mncg, int numbuckets)
{
    int i, j;
    float step=(mxcg-mncg)/(float)numbuckets;
    float *bucketlimarr=malloc((numbuckets-1)*sizeof(float));
    int *bucketarr=calloc(numbuckets, sizeof(int));
    bucketlimarr[0]=step+mncg;
    for(i=1;i<numbuckets-1;++i) 
        bucketlimarr[i]=bucketlimarr[i-1]+step;

    for(i=0;i<sz;++i)
        if(sqisz[i].cgp>=bucketlimarr[numbuckets-2]) {
            bucketarr[numbuckets-1]++;
            continue;
        } else {
            for(j=0;j<numbuckets-1;++j)
                if(sqisz[i].cgp<bucketlimarr[j]) {
                    bucketarr[j]++;
                    break;
                }
        }
    free(bucketlimarr);
    return bucketarr;
}

int *hist_tsz(i_s *sqisz, int sz, size_t mxtsz, size_t mntsz, int numbuckets)
{
    int i, j;
    float step=(float)(mxtsz-mntsz)/numbuckets;
    float *bucketlimarr=malloc((numbuckets-1)*sizeof(float));
    int *bucketarr=calloc(numbuckets, sizeof(int));
    bucketlimarr[0]=step+(float)mntsz;
    for(i=1;i<numbuckets-1;++i) 
        bucketlimarr[i]=bucketlimarr[i-1]+step;

    for(i=0;i<sz;++i)
        if(sqisz[i].tsz>=bucketlimarr[numbuckets-2]) {
            bucketarr[numbuckets-1]++;
            continue;
        } else {
            for(j=0;j<numbuckets-1;++j) /* yes, the last bucketlim is being reused! this time, to find values below it! */
                if(sqisz[i].tsz<bucketlimarr[j]) {
                    bucketarr[j]++;
                    break;
                }
        }
    free(bucketlimarr);
    return bucketarr;
}

int *hist_ambs(i_s *sqisz, int sz, unsigned mxamb, unsigned mnamb, int numbuckets)
{
    int i, j;
    float step=(float)(mxamb-mnamb)/numbuckets;
    float *bucketlimarr=malloc((numbuckets-1)*sizeof(float));
    int *bucketarr=calloc(numbuckets, sizeof(int));
    bucketlimarr[0]=step+(float)mnamb;
    for(i=1;i<numbuckets-1;++i) 
        bucketlimarr[i]=bucketlimarr[i-1]+step;

    for(i=0;i<sz;++i)
        if(sqisz[i].ambano[0]>=bucketlimarr[numbuckets-2]) {
            bucketarr[numbuckets-1]++;
            continue;
        } else {
            for(j=0;j<numbuckets-2;++j)
                if(sqisz[i].ambano[0]<bucketlimarr[j]) {
                    bucketarr[j]++;
                    break;
                }
        }
    free(bucketlimarr);
    return bucketarr;
}

void prti_s(i_s *sqisz, int sz, float *mxcg, float *mncg)
{
    int i;
    char *sqgood;
    *mxcg=.0;
    *mncg=1.;

    size_t tsz;
    for(i=0;i<sz;++i) {
        if(sqisz[i].ambano[1] != 0)
            sqgood="AnoSQ";
        else
            sqgood="SQ";
        tsz = sqisz[i].sy[0] + sqisz[i].sy[1];
        sqisz[i].cgp=(float)sqisz[i].sy[0]/tsz;
        if(sqisz[i].cgp>*mxcg)
            *mxcg=sqisz[i].cgp;
        if(sqisz[i].cgp<*mncg)
            *mncg=sqisz[i].cgp;

        printf("| %s#%i=TOT:%zu CG:%.4f ", sqgood, i, sqisz[i].tsz, sqisz[i].cgp);
    }
    printf("|\n"); 
}

int main(int argc, char *argv[])
{
    /* argument accounting: remember argc, the number of arguments, _includes_ the executable */
    if(argc!=2) {
        printf("Error. Pls supply 1 argument:; fasta file name.\n");
        exit(EXIT_FAILURE);
    }
    FILE *fin;
    if(!(fin=fopen(argv[1], "r")) ) { /*s houdl one check the extension of the fasta file ? */
        printf("Error. Cannot open the presented filename\n");
        exit(EXIT_FAILURE);
    }

    char IGLINE=0, begline=1;
    size_t lidx=0, mxtsz=0, mntsz=0XFFFFFFFFFFFFFFFF;
    unsigned mxamb=0, mnamb=0xFFFFFFFF;

    int i, c, sqidx=-1; /* this is slightly dangerous, you need very much to knwo what you're doing */
    int gbuf=GBUF;
    i_s *sqisz=malloc(gbuf*sizeof(i_s));
    sqisz.intsqbuf=SBUF;
    intsq.intsq=calloc(SBUF, sizeof(char));
    char whatint;

    while( ( (c = fgetc(fin)) != EOF) ) {
        if(c =='\n') {
            IGLINE=0;
            begline=1;
            lidx++;
        } else if( (begline==1) & (c == '>') ) { /* this condition catche sthe beginning of a new sequence, and uses it to prepare the nextsequence.*/
            IGLINE =1;
            begline=0; 
            if(sqidx>=0) { /* multifasta detected ... we will only consider the first onee */
                if(sqisz[sqidx].tsz > mxtsz)
                    mxtsz = sqisz[sqidx].tsz;
                if(sqisz[sqidx].tsz < mntsz)
                    mntsz = sqisz[sqidx].tsz;
                if(sqisz[sqidx].ambano[0] > mxamb)
                    mxamb = sqisz[sqidx].ambano[0];
                if(sqisz[sqidx].ambano[0] < mnamb)
                    mnamb = sqisz[sqidx].ambano[0];
            }

            sqidx++;
            if(sqidx==gbuf) {
                gbuf+=GBUF;
                sqisz=realloc(sqisz, gbuf*sizeof(i_s));
            }
            sqisz[sqidx].idx=sqidx;
            sqisz[sqidx].tsz=0;
            for(i=0;i<SSZ;++i)
                sqisz[sqidx].sy[i]=0;
            for(i=0;i<2;++i)
                sqisz[sqidx].ambano[i]=0;
        } else if (IGLINE==0) {
            sqisz[sqidx].tsz++;
            if(sqidx[sqidx].tsz==sqidx[sqidx].intsqbuf) {
                sqidx[sqidx].intsqbuf+=SBUF;
                sqisz[sqidx].intsq=realloc(sqisz[sqidx].intsq, sqidx[sqidx].intsqbuf*sizeof(char));
            }
            switch(c) {
                case 'A': case 'a':
                    whatint=1; break;
                case 'C': case 'c':
                    whatint=2; break;
                case 'G': case 'g':
                    whatint=3; break;
                case 'T': case 't':
                    whatint=4; break;
                case 'R': case 'r':
                    whatint=5; break;
                case 'Y': case 'y':
                    whatint=6; break;
                case 'K': case 'k': /* the ketos */
                    whatint=7; break;
                case 'M': case 'm': /* the aminoids */
                    whatint=8; break;
                case 'S': case 's':
                    whatint=9; break;
                case 'W': case 'w':
                    whatint=10; break;
                case 'B': case 'b':
                    whatint=11; break;
                case 'D': case 'd':
                    whatint=12; break;
                case 'H': case 'h':
                    whatint=13; break;
                case 'V': case 'v':
                    whatint=14; break;
                case 'N': case 'n':
                    whatint=15; break;
                case '-':
                    whatint=16; break;
                default:
                    whatint=17; /* unknown this means your fasta file is naff. */
                    sqisz[sqidx].ambano[1]++; break;
            }
        }

    }
    fclose(fin);
    /* postprocessing on the final sequence */
    if(sqisz[sqidx].tsz > mxtsz)
        mxtsz = sqisz[sqidx].tsz;
    if(sqisz[sqidx].tsz < mntsz)
        mntsz = sqisz[sqidx].tsz;
    if(sqisz[sqidx].ambano[0] > mxamb)
        mxamb = sqisz[sqidx].ambano[0];
    if(sqisz[sqidx].ambano[0] < mnamb)
        mnamb = sqisz[sqidx].ambano[0];

    unsigned numsq=sqidx+1, numano=0;
    for(i=0;i<numsq;++i) 
        if(sqisz[i].ambano[1])
            numano++;
    sqisz=realloc(sqisz, numsq*sizeof(i_s));
    float mxcg, mncg;
    prti_s(sqisz, numsq, &mxcg, &mncg);
    /* the summary comes at the end because otherwise, with many sequences, it goes off-screen */
    printf("Number of sequences: %i, Mxsz= %zu, Minsz= %zu, MaxCG=%.4f MinCG=%.4f Mxamb=%u Mnamb=%u. #AnoSQ=%i\n", numsq, mxtsz, mntsz, mxcg, mncg, mxamb, mnamb, numano);

    /* OK tsz histo first */
    int numbuckets=HISTBUCKETSZ;
    int *histosz=hist_tsz(sqisz, numsq, mxtsz, mntsz, numbuckets);
    prthist("Seqsz", histosz, numbuckets);
    int *histocg=hist_cg(sqisz, numsq, mxcg, mncg, numbuckets);
    prthist("cgpart", histocg, numbuckets);

    free(histosz);
    free(histocg);
    free(sqisz);
    return 0;
}
