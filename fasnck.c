/* nfasack.c DNA Fasta file sanity check: principally it will say 
 * if it a multisequence file */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GBUF 4
#define SSZ 2 /* CG count, first, AT count second, third are the anomalous characters */
#define HISTBUCKETSZ 10

typedef struct /* i_s; sequence index and number of symbols */
{
    unsigned int idx;
    size_t sylen; /* this is the precise symbol length of the sequence */
    size_t sy[SSZ]; /* used to hold counts of symbols */
    float cgp;
    unsigned ambano[2]; /* number of ambiguous symbots (first), then number of anomalous symbols */
} i_s; /* sequence index and number of symbols */

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

int *hist_sylen(i_s *sqisz, int sz, size_t mxsylen, size_t mnsylen, int numbuckets)
{
    int i, j;
    float step=(float)(mxsylen-mnsylen)/numbuckets;
    float *bucketlimarr=malloc((numbuckets-1)*sizeof(float));
    int *bucketarr=calloc(numbuckets, sizeof(int));
    bucketlimarr[0]=step+(float)mnsylen;
    for(i=1;i<numbuckets-1;++i) 
        bucketlimarr[i]=bucketlimarr[i-1]+step;

    for(i=0;i<sz;++i)
        if(sqisz[i].sylen>=bucketlimarr[numbuckets-2]) {
            bucketarr[numbuckets-1]++;
            continue;
        } else {
            for(j=0;j<numbuckets-1;++j) /* yes, the last bucketlim is being reused! this time, to find values below it! */
                if(sqisz[i].sylen<bucketlimarr[j]) {
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

void prthist(char *histname, int *bucketarr, int numbuckets, unsigned numsq, size_t mxsylen, size_t mnsylen)
{
    int i;
    printf("Sqlen %d-bin hstgrm for: %-24.24s (totsqs=%04u): ", numbuckets, histname, numsq); 
    printf("minlen=%3zu<-", mnsylen); 
    for(i=0;i<numbuckets;++i) 
        printf("| %i ", bucketarr[i]);
    printf("|->maxlen=%zu\n", mxsylen); 
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

        printf("| %s#%i=TOT:%zu CG:%.4f ", sqgood, i, sqisz[i].sylen, sqisz[i].cgp);
    }
    printf("|\n"); 
}

void la_prti_s(i_s *sqisz, int sz, float *mxcg, float *mncg, char *titlestr) /* prints onto a beamer table, also calculates mx amnd min cg in passsing */
{
    int i, cols=4;
    printf("\\begin{frame}\n\t\\frametitle{%s}\n", titlestr);
    char *h0[4]= {"Seq Idx", "Presence AmbSymbs", "Seq Length", "\\% CG Content"};
    size_t tsz;
    char *sqgood;

    printf("\\begin{center}\n");
    printf("\\small\n");
    printf("    \\begin{tabular}{| l | r | r | r |}\n");
    printf("        \\hline\n               ");
    for(i=0;i<cols;++i) {
        printf("%s", h0[i]);
        (i != cols-1)? printf(" & ") : printf(" \\\\\n");
    }
    printf("        \\hline\n");

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

        printf("%i & %s & %zu & %.3f \\\\\n", i, sqgood, sqisz[i].sylen, sqisz[i].cgp);
    }

    printf("        \\hline\n");
    printf("    \\end{tabular}\n");
    printf("\\end{center}\n");

    printf("\\end{frame}\n\n");

    return;
}

void la_prti_s2(i_s *sqisz, int sz, float *mxcg, float *mncg, char *titlestr) /* version2 moidifies printing onto a beamer table, to generic STDOUT table */
{
    int i, cols=4;
    char *h0[4]= {"Seq Idx", "Presence AmbSymbs", "Seq Length", "\\% CG Content"};
    size_t tsz;
    char *sqgood;

    printf("--------------------------------\n");
    for(i=0;i<cols;++i) {
        printf("%s", h0[i]);
        (i != cols-1)? printf(" & ") : printf(" \\\\\n");
    }
    printf("--------------------------------\n");

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

        printf("%i | %s | %zu | %.3f \n", i, sqgood, sqisz[i].sylen, sqisz[i].cgp);
    }

    printf("--------------------------------\n");

    return;
}

void tikz_prti_s(i_s *sqisz, int sz, int *histosz, int numbuckets, char *titlestr) /* prints onto a beamer table, also calculates mx amnd min cg in passsing */
{
    int j;
    printf("\\begin{frame}\n\t\\frametitle{%s}\n", titlestr);
    char *h[2]= {"Seq Len", "Seq Quan"};

    printf("\\begin{tikzpicture}[scale=1.0]\n");
    printf("    \\begin{axis}[\n");

    printf("     xlabel=%s,\n", h[0]);
    printf("     ylabel=%s]\n", h[1]);

    printf("    \\addplot[color=%s] coordinates {\n", "blue");
    for(j=0;j<numbuckets;++j)
        printf("       (%i,%i)\n", j, histosz[j]);
    printf("    };\n");

    printf("    \\end{axis}\n");
    printf("\\end{tikzpicture}\n");
    printf("\\end{frame}\n\n");

    return;
}

int main(int argc, char *argv[])
{
    /* argument accounting: remember argc, the number of arguments, _includes_ the executable */
    if(argc==1) {
        printf("Error. Pls supply 1+ arguments: A sequence of multi-fasta filenames. \n");
        exit(EXIT_FAILURE);
    }
    /* general declarations */
    FILE *fin;
    char IGLINE, begline;
    size_t lidx, mxsylen, mnsylen;
    unsigned mxamb, mnamb;
    int i, j, c, sqidx;
    int gbuf;
    i_s *sqisz;
    int whatint;
    int numbuckets;
    int *histosz;
    unsigned numsq, numano;

    for(j=1;j<argc;++j) {

        if(!(fin=fopen(argv[j], "r")) ) { /*s houdl one check the extension of the fasta file ? */
            printf("Error. Cannot open \"%s\" file.\n", argv[j]);
            exit(EXIT_FAILURE);
        }

        IGLINE=0, begline=1;
        lidx=0, mxsylen=0, mnsylen=0XFFFFFFFFFFFFFFFF;
        mxamb=0, mnamb=0xFFFFFFFF;

        sqidx=-1; /* this is slightly dangerous, you need very much to knwo what you're doing */
        gbuf=GBUF;
        sqisz=malloc(gbuf*sizeof(i_s));
        for(i=gbuf-GBUF;i<gbuf;++i) {
            sqisz[i].ambano[0]=0;
            sqisz[i].ambano[1]=0;
        }
        whatint=0; /* needs explanation */

        while( ( (c = fgetc(fin)) != EOF) ) {
            if(c =='\n') {
                IGLINE=0;
                begline=1;
                lidx++;
            } else if( (begline==1) & (c == '>') ) { /* this condition catches the beginning of a new sequence, and uses it to prepare the nextsequence.*/
                IGLINE =1;
                begline=0; 
                if(sqidx>=0) { /* chancing my arm here ... operating on the past sequence */
                    if(sqisz[sqidx].sylen > mxsylen)
                        mxsylen = sqisz[sqidx].sylen;
                    if(sqisz[sqidx].sylen < mnsylen)
                        mnsylen = sqisz[sqidx].sylen;
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
                sqisz[sqidx].sylen=0;
                for(i=0;i<SSZ;++i)
                    sqisz[sqidx].sy[i]=0;
                for(i=0;i<2;++i)
                    sqisz[sqidx].ambano[i]=0;
            } else if (IGLINE==0) {
                sqisz[sqidx].sylen++;
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
                }
            }
            if( (whatint == 2) || (whatint == 3) ) {
                sqisz[sqidx].sy[0]++;
                sqisz[sqidx].ambano[1]++;
            } else if (whatint < 5) {
                sqisz[sqidx].sy[1]++;
                sqisz[sqidx].ambano[1]++;
            } else 
                sqisz[sqidx].ambano[0]++;
        }
        fclose(fin);
        /* postprocessing on the final sequence */
        if(sqisz[sqidx].sylen > mxsylen)
            mxsylen = sqisz[sqidx].sylen;
        if(sqisz[sqidx].sylen < mnsylen)
            mnsylen = sqisz[sqidx].sylen;
        if(sqisz[sqidx].ambano[0] > mxamb)
            mxamb = sqisz[sqidx].ambano[0];
        if(sqisz[sqidx].ambano[0] < mnamb)
            mnamb = sqisz[sqidx].ambano[0];

        numsq=sqidx+1, numano=0;
        for(i=0;i<numsq;++i) {
            if(sqisz[i].ambano[1])
                numano++;
        }
        sqisz=realloc(sqisz, numsq*sizeof(i_s));
        float mxcg, mncg;
        la_prti_s2(sqisz, numsq, &mxcg, &mncg, "Table of Sequence Lengths");
        //  prti_s(sqisz, numsq, &mxcg, &mncg);

        /* OK sylen histo first */
        numbuckets=HISTBUCKETSZ;
        histosz=hist_sylen(sqisz, numsq, mxsylen, mnsylen, numbuckets);
        // tikz_prti_s(sqisz, numsq, histosz, numbuckets, "Table of Sequence Lengths");
        prthist(argv[j], histosz, numbuckets, numsq, mxsylen, mnsylen);
        //    int *histocg=hist_cg(sqisz, numsq, mxcg, mncg, numbuckets);
        // prthist("cgpart", histocg, numbuckets);
        /* the summary comes at the end because otherwise, with many sequences, it goes off-screen */
        fprintf(stderr, "Number of sequences: %i, Mxsz= %zu, Minsz= %zu, MaxCG=%.4f MinCG=%.4f Mxamb=%u Mnamb=%u. #AnoSQ=%i\n", numsq, mxsylen, mnsylen, mxcg, mncg, mxamb, mnamb, numano);

        free(histosz);
        //     free(histocg);
        free(sqisz);
    }
    return 0;
}
