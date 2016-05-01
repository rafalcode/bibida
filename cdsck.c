/* nfasack.c DNA Fasta file sanity check: principally it will say 
 * if it a multisequence file */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef DBG
#define GBUF 4
#else
#define GBUF 128
#endif
#define SSZ 2 /* CG count, first, AT count second, third are the anomalous characters */
#define HISTBUCKETSZ 10

#define CONDREALLOC(x, b, c, a, t); \
    if((x)==((b)-1)) { \
        (b) += (c); \
        (a)=realloc((a), (b)*sizeof(t)); \
        memset((a)+(b)-(c), '\0', (c)*sizeof(t)); \
    }

const char *stcod="ATG";
const char *ecod1="TAG";
const char *ecod2="TAA";
const char *ecod3="TGA";

typedef struct /* onefa */
{
    char *id;
    char *sq;
    unsigned idz, sqz;
} onefa;

typedef struct /* i_s; sequence index and number of symbols */
{
    unsigned int idx;
    size_t sylen; /* this is the precise symbol length of the sequence */
    size_t sy[SSZ]; /* used to hold counts of symbols */
    float cgp;
    unsigned ambano[2]; /* number of ambiguous symbols (first), then number of anomalous symbols */
} i_s; /* sequence index and number of symbols */

void prtfa(onefa *fac)
{
    printf(">");
    printf("%s\n", fac->id);
    printf("%s\n", fac->sq);
}

void prtfaf(onefa *fac, FILE *fp)
{
    fprintf(fp, ">");
    fprintf(fp, "%s\n", fac->id);
    fprintf(fp, "%s\n", fac->sq);
}

void prtfa2(onefa *fac)
{
    int i;
    printf("SQZ=%d:", fac->sqz);
    for(i=0;i<3;++i) 
        putchar(fac->sq[i]);
    printf("\n"); 
}

void ck(onefa *fac)
{
    char repstr[128]={0};
    if(((fac->sqz-1)%3) != 0)
        strcat(repstr, "E_%3 ");
    if( (strncmp(fac->sq, stcod, 3*sizeof(char))) )
        strcat(repstr, "E_stcod ");
    char *lcod=fac->sq+(fac->sqz-4);
    if( (strncmp(lcod, ecod1, 3*sizeof(char))) & (strncmp(lcod, ecod2, 3*sizeof(char))) & (strncmp(lcod, ecod3, 3*sizeof(char))) )
        strcat(repstr, "E_ecod ");

    if(repstr[0])
        printf("%s: %s\n", fac->id, repstr); 
}

int cki(onefa *fac)
{
    int i=0;
    if(((fac->sqz-1)%3) != 0)
        i++;
    if( (strncmp(fac->sq, stcod, 3*sizeof(char))) )
        i++;
    char *lcod=fac->sq+(fac->sqz-4);
    if( (strncmp(lcod, ecod1, 3*sizeof(char))) & (strncmp(lcod, ecod2, 3*sizeof(char))) & (strncmp(lcod, ecod3, 3*sizeof(char))) )
        i++;
    return i;
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

int main(int argc, char *argv[])
{
    /* argument accounting: remember argc, the number of arguments, _includes_ the executable */
    if(argc==1) {
        printf("Error. Pls supply 1+ arguments: A sequence of multi-fasta filenames. \n");
        exit(EXIT_FAILURE);
    }
    /* general declarations */
    FILE *fin, *fout;
    char *tp, *foutname=calloc(128, sizeof(char));
    char insertstr[10]="_filtered";
    char IGLINE, begline;
    size_t lidx, mxsylen, mnsylen;
    unsigned mxamb, mnamb;
    int i, j, c, sqidx;
    int gbuf;
    i_s *sqisz;
    int whatint;
    unsigned numsq, numano;
    onefa fac;
    unsigned ibf=GBUF, sbf=GBUF;
    fac.id=calloc(ibf, sizeof(char));
    fac.sq=calloc(sbf, sizeof(char));
    size_t fnsz;
    int ididx=0;
    int retcki=0;
    unsigned numfiltered=0;


    for(j=1;j<argc;++j) {

        if(!(fin=fopen(argv[j], "r")) ) { /*should one check the extension of the fasta file ? */
            printf("Error. Cannot open \"%s\" file.\n", argv[j]);
            exit(EXIT_FAILURE);
        }
        
        fnsz=1+strlen(argv[j]);
        foutname=realloc(foutname, (fnsz+10)*sizeof(char));
        tp=strrchr(argv[j], '.');
        sprintf(foutname, "%.*s%s%s", (int)(tp-argv[j]), argv[j], insertstr, tp);
        fout=fopen(foutname, "w");
        IGLINE=0, begline=1;
        lidx=0, mxsylen=0, mnsylen=0XFFFFFFFFFFFFFFFF;
        mxamb=0, mnamb=0xFFFFFFFF;

        sqidx=-1; /* this is slightly dangerous, you need very much to know what you're doing */
        gbuf=GBUF;
        sqisz=malloc(gbuf*sizeof(i_s));
        for(i=gbuf-GBUF;i<gbuf;++i) {
            sqisz[i].ambano[0]=0;
            sqisz[i].ambano[1]=0;
        }
        whatint=0; /* needs explanation */
        ididx=0;
        retcki=0;
        numfiltered=0;

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

                    CONDREALLOC(ididx, ibf, GBUF, fac.id, char);
                    fac.id[ididx]='\0';
                    CONDREALLOC(sqisz[sqidx].sylen, sbf, GBUF, fac.sq, char);
                    fac.sq[sqisz[sqidx].sylen]='\0';
                    fac.idz=1+ididx;
                    fac.sqz=1+sqisz[sqidx].sylen;
                    // prtfa2(&fac);
                    retcki=cki(&fac);
                    if(!retcki) {
                        prtfaf(&fac, fout);
                        numfiltered++;
                    }
                }

                sqidx++;
                if(sqidx==gbuf) {
                    gbuf+=GBUF;
                    sqisz=realloc(sqisz, gbuf*sizeof(i_s));
                }
                sqisz[sqidx].idx=sqidx;

                /* resetting stuff */
                sqisz[sqidx].sylen=0;
                ididx=0;
                ibf=GBUF;
                sbf=GBUF;
                fac.id=realloc(fac.id, ibf*sizeof(char));
                memset(fac.id, '\0', ibf*sizeof(char));
                fac.sq=realloc(fac.sq, sbf*sizeof(char));
                memset(fac.sq, '\0', sbf*sizeof(char));
                for(i=0;i<SSZ;++i)
                    sqisz[sqidx].sy[i]=0;
                for(i=0;i<2;++i)
                    sqisz[sqidx].ambano[i]=0;
            } else if (IGLINE==1) {
                CONDREALLOC(ididx, ibf, GBUF, fac.id, char);
                fac.id[ididx++]=c;
            } else if (IGLINE==0) {
                CONDREALLOC(sqisz[sqidx].sylen, sbf, GBUF, fac.sq, char);
                fac.sq[sqisz[sqidx].sylen]=c;
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

        /* the last sequence */
        CONDREALLOC(ididx, ibf, GBUF, fac.id, char);
        fac.id[ididx]='\0';
        CONDREALLOC(sqisz[sqidx].sylen, sbf, GBUF, fac.sq, char);
        fac.sq[sqisz[sqidx].sylen]='\0';
        fac.idz=1+ididx;
        fac.sqz=1+sqisz[sqidx].sylen;
       // prtfa2(&fac);
        retcki=cki(&fac);
        if(!retcki) {
            numfiltered++;
            prtfaf(&fac, fout);
        }
        fclose(fout);
        
        /* free */
        free(fac.id);
        fac.id=NULL;
        free(fac.sq);
        fac.sq=NULL;

        numsq=sqidx+1, numano=0;
        for(i=0;i<numsq;++i) {
            if(sqisz[i].ambano[1])
                numano++;
        }
        sqisz=realloc(sqisz, numsq*sizeof(i_s));
        printf("%s/numseq=%u/umfiltered=%u/%%filt=%.1f\n", argv[j], numsq, numfiltered, 100*((float)numfiltered/numsq)); 
        // prti_s(sqisz, numsq, &mxcg, &mncg);
        /* the summary comes at the end because otherwise, with many sequences, it goes off-screen */
        // fprintf(stderr, "Number of sequences: %i, Mxsz= %zu, Minsz= %zu, MaxCG=%.4f MinCG=%.4f Mxamb=%u Mnamb=%u. #AnoSQ=%i\n", numsq, mxsylen, mnsylen, mxcg, mncg, mxamb, mnamb, numano);

        free(sqisz);
    }
    free(foutname);
    return 0;
}
