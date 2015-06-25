/* fastitch.c stitches up fasta files */
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#ifdef DBG
#define GBUF 6
#define PDCHAR ' ' /* put this up at top, what the padding character will be */
#define MINPADNLEN 3 /* normally we also want to apply a minimum N-pad so that the short reads can map independently to the stitched contig components */
#else
#define GBUF 128
#define PDCHAR 'N' /* put this up at top, what the padding character will be */
#define MINPADNLEN 300 /* normally we also want to apply a minimum N-pad so that the short reads can map independently to the stitched contig components */
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

#define CONDREALLOC2(x, b, c, a, t); \
    if((x)==(b)) { \
        (b) += (c); \
        (a)=realloc((a), (b)*sizeof(t)); \
    }

#define SIACONDREALLOC(x, b, c, a, i); \
    if((x)>=((b)-2)) { \
        (b) += (c); \
        (a)=realloc((a), (b)*sizeof(sia_t)); \
        for((i)=((b)-(c));(i)<(b);(++(i))) { \
            (a[i].sisz)=0; \
            (a[i].tssz)=0; \
            (a[i].mxssz)=0; \
            (a[i].sibf)=(b); \
            (a[i].si)=malloc((a[i].sibf)*sizeof(int)); \
        } \
    }


typedef struct /* uoa Unique Occurence Array type: the UO here refers to the length of the sequence: i.e. the occurence of a certain sequnce length. You'll need to remember that :-) . */
{
    unsigned uo; /* the unique occurence this array entry refers to */
    int *uoids; /* the ids of the sequences for this unqiue occurence */
    int uoisz; /* the actual size of the array */
} uoa_t;

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

typedef struct /* i_sa type */
{
    i_s *is;
    int numsq;
    unsigned totbases;
    unsigned mx, mn;
} i_sa;

typedef struct /* sia_t */
{
    int *si; /* sequence idx */
    unsigned sibf;
    unsigned sisz;
    unsigned tssz; /* actual accum size of the sequences */
    unsigned mxssz; /* the max uo, i.e. max size of seq compoenents */
} sia_t;

i_sa *crea_i_sa(void)
{
    i_sa *sqia=malloc(1*sizeof(i_sa));
    sqia->is=NULL;
    sqia->numsq=0;
    sqia->totbases=0;
    // for(i=0;i<sqia->numsq;i++) // thinking of it
    sqia->mn=0xFFFFFFFF;
    sqia->mx=0;
    return sqia;
}

void free_i_sa(i_sa **sqia_)
{
    int i;
    i_sa *sqia=*sqia_;
    for(i=0; i<sqia->numsq;++i) {
        free(sqia->is[i].id);
        free(sqia->is[i].sq);
    }
    free(sqia->is);
    free(sqia);
    sqia=NULL;
    sqia_=NULL;
}

void usage(char *executable)
{
    printf("%s, a program to stitch a fragmented fasta file (usually a reference file).\n", executable);
    printf("Usage: %s <one argument: the (usually) reference fasta filename.\n", executable);
    printf("Explanation:\n");
    printf("\tThis is an improved second version (i.e. beyond the hack that the first one was)\n");
    printf("\tthat groups the sequences which all have the same size, and merges them together in a certain manner.\n");
    printf("\tUsing N's of course. There's quite a little science involved in uniformizing lengths in this way, and\n");
    printf("\tthere's quite a little way to go. But, if you are using it, best to run \"fasnck\" on the stitched result\n");
    printf("\tto see if you are getting anywhere fast.\n");
}

int cmpuoabyo(const void *a, const void *b) /* compare uoa by occurence */
{
    uoa_t *ua = (uoa_t*)a; /* cast our void! */
    uoa_t *ub = (uoa_t*)b; /* cast our void! */
    return ua->uo  - ub->uo; /* integer comparison: returns positive if b > a and nagetive if a > b: i.e. highest values first */
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

i_s *delthese(i_s **sqi_, unsigned *nsq, unsigned *idstodel, unsigned idstodeln) /* returns a new i_s with the ids that are deleted from the main */
{
    int i;
    int numsq=*nsq;
    i_s *sqi=*sqi_;
    i_s *sqi2=calloc(idstodeln, sizeof(i_s));

    for(i=0; i<idstodeln; ++i) {
        sqi2[i].sylen = sqi[idstodel[i]].sylen;
        sqi2[i].sqz = sqi[idstodel[i]].sqz;
        sqi2[i].sq=calloc(sqi2[i].sqz, sizeof(char));
        memcpy(sqi2[i].sq, sqi[idstodel[i]].sq, sqi[2].sqz*sizeof(char));
        sqi2[i].idz = sqi[idstodel[i]].idz;
        sqi2[i].id=calloc(sqi2[i].idz, sizeof(char));
        memcpy(sqi2[i].id, sqi[idstodel[i]].id, sqi[2].idz*sizeof(char));
    }
    return sqi2;
    *nsq=numsq;
}

void prtasf(i_s *sqi, int sz, FILE *fp) /* print all sequences to file */
{
    int i;
    for(i=0;i<sz;++i) {
        fprintf(fp, ">%s\n", sqi[i].id);
        fprintf(fp, "%s\n", sqi[i].sq);
    }
}

void prtseq0(i_s *sqi, int sz, FILE *fp, int *whichids, int whichidsz) /* straight print of a selected array of sequences */
{
    int i, j;
    for(i=0;i<whichidsz;++i) {
        fprintf(fp, ">");
        for(j=0;j<sqi[whichids[i]].idz-1;j++)
            fprintf(fp, "%c", sqi[whichids[i]].id[j]);
        fprintf(fp, "\n"); 
        for(j=0;j<sqi[whichids[i]].sqz-1;j++)
            fprintf(fp, "%c", sqi[whichids[i]].sq[j]);
        fprintf(fp, "\n"); 
    }
}

void prtseq1(i_s *sqi, int sz, FILE *fp, unsigned uo, int *whichids, int whichidsz) /* print several sequences into one, padded by Ns */
{
    int i, j;
    char paddingchar=PDCHAR;
    int minnlen=MINPADNLEN;
    unsigned mxclen=uo;
    int npadlen = (mxclen > minnlen)? mxclen : minnlen;

    for(i=0;i<whichidsz;++i) {
        if(i==0) { /* the idline used will be that of the first sequence */
            fprintf(fp, ">");
            for(j=0;j<sqi[whichids[i]].idz-1;j++)
                fprintf(fp, "%c", sqi[whichids[i]].id[j]);
            fprintf(fp, "\n"); 
        }
        for(j=0;j<sqi[whichids[i]].sqz-1;j++)
            fprintf(fp, "%c", sqi[whichids[i]].sq[j]);
        if( i != whichidsz-1) 
            for(j=0;j<npadlen;++j) 
                fputc(paddingchar, fp);
    }
    fprintf(fp, "\n"); 
}


void prtsiaele(i_sa *sqia, FILE *fp, int *whichuois, int whichuoisz, unsigned mxclen) /* print a sia element ... actually sia members used and they refer directly to sqia->si */
{
    int i, k;
    char paddingchar=PDCHAR;
    int minnlen=MINPADNLEN;
    int npadlen = (mxclen > minnlen)? mxclen : minnlen;

    /* really should be modifying the ID string of the first sequence */
    for(i=0;i<whichuoisz;++i) { /* all these have to go in one sequence */
        if(i==0) {
            fprintf(fp, ">");
            for(k=0;k<sqia->is[whichuois[i]].idz-1;k++)
                fprintf(fp, "%c", sqia->is[whichuois[i]].id[k]);
            fprintf(fp, "%s", (whichuoisz != 1)? "[NB:MGDSEQ]\n" : "\n"); 
        }
        for(k=0;k<sqia->is[whichuois[i]].sqz-1;k++)
            fprintf(fp, "%c", sqia->is[whichuois[i]].sq[k]);
        if(i != whichuoisz-1)
            for(k=0;k<npadlen;++k) 
                fputc(paddingchar, fp);
    }
    fprintf(fp, "\n"); 
}

uoa_t *uniquelens(i_sa *sqia, int *uoasz_)
{
    unsigned char new;
    unsigned i, j;
    unsigned uoabuf=GBUF;
    int uoasz=0;
    uoa_t *uoa=calloc(uoabuf, sizeof(uoa_t));
    for(i=0; i<sqia->numsq;++i) {
        new=1;
        for(j=0; j<uoasz;++j) {
            if(uoa[j].uo == sqia->is[i].sylen) {
                uoa[j].uoisz++;
                uoa[j].uoids=realloc(uoa[j].uoids, uoa[j].uoisz*sizeof(unsigned));
                uoa[j].uoids[uoa[j].uoisz-1] = i;
                // uoa[uoasz-1].uoids[uoa[j].uoisz-1] = sqia->is[i].idx;
                new=0;
                break;
            }
        }
        if(new) {
            uoasz++;
            CONDREALLOC2(uoasz, uoabuf, GBUF, uoa, uoa_t);
            uoa[uoasz-1].uo = sqia->is[i].sylen;
            uoa[uoasz-1].uoisz = 1;
            uoa[uoasz-1].uoids=NULL;
            uoa[uoasz-1].uoids=realloc(uoa[j].uoids, uoa[j].uoisz*sizeof(unsigned));
            uoa[uoasz-1].uoids[uoa[j].uoisz-1] = i;
        }
    }

    qsort(uoa, uoasz, sizeof(uoa_t), cmpuoabyo);
#ifdef DBG
    printf("number of different sequence lengths: %i\n", uoasz);
    for(j=0; j<uoasz;++j) {
        printf("%i (%u): ", uoa[j].uo, uoa[j].uoisz);
        for(i=0;i<uoa[j].uoisz;++i) 
            printf("%i ", uoa[j].uoids[i]);
        printf("\n"); 
    }
#endif
    *uoasz_=uoasz;
    return uoa;
}

void fprtuoahist(char *fname, uoa_t *uoa, int uoasz)
{
    int j;
    FILE *fhist=fopen(fname, "w");
    for(j=0; j<uoasz;++j)
        fprintf(fhist, "%i %u\n", uoa[j].uo, uoa[j].uoisz);
    fclose(fhist);
}

i_sa *faf_to_i_s(char *fafname) /* fasta file to i_s data structure */
{
    /* general declarations */
    FILE *fin;
    char IGLINE, begline;
    size_t lidx;
    int i, c, sqidx;
    int gbuf;
    int cx=0; /* character index */
    i_sa *sqia=crea_i_sa();

    if(!(fin=fopen(fafname, "r")) ) { /* should one check the extension of the fasta file? */
        printf("Error. Cannot open \"%s\" file.\n", fafname);
        exit(EXIT_FAILURE);
    }

    IGLINE=0, begline=1, lidx=0, sqidx=-1; /* that last one is slightly dangerous, you need to know what you're doing */
    gbuf=GBUF;
    sqia->is=malloc(gbuf*sizeof(i_s));
    for(i=gbuf-GBUF;i<gbuf;++i) {
        sqia->is[i].ibf=GBUF;
        sqia->is[i].sbf=GBUF;
        sqia->is[i].id=calloc(sqia->is[i].ibf, sizeof(char));
        sqia->is[i].sq=calloc(sqia->is[i].sbf, sizeof(char));
    }
    cx=0;

    while( ( (c = fgetc(fin)) != EOF) ) {
        if(c =='\n') {
            IGLINE=0;
            begline=1;
            lidx++;
        } else if( (begline==1) & (c == '>') ) { /* this condition catches the beginning of a new sequence, and uses it to prepare the nextsequence.*/
            IGLINE =1;
            begline=0; 
            if(sqidx>=0) { /* we're not interested in the first run, when sqidx=-1. Our work on the sqi array is "retrospective" */

                CONDREALLOC_(cx, sqia->is[sqidx].ibf, GBUF, sqia->is[sqidx].id, char);
                sqia->is[sqidx].id[cx]='\0';
                CONDREALLOC(sqia->is[sqidx].sylen, sqia->is[sqidx].sbf, GBUF, sqia->is[sqidx].sq, char);
                sqia->is[sqidx].sq[sqia->is[sqidx].sylen]='\0';
                sqia->is[sqidx].idz=1+cx;
                sqia->is[sqidx].sqz=1+sqia->is[sqidx].sylen;
                sqia->totbases += sqia->is[sqidx].sylen;
                if(sqia->is[sqidx].sylen > sqia->mx)
                    sqia->mx = sqia->is[sqidx].sylen;
                if(sqia->is[sqidx].sylen < sqia->mn)
                    sqia->mn = sqia->is[sqidx].sylen;
            }

            sqidx++;
            if(sqidx==gbuf) {
                gbuf+=GBUF;
                sqia->is=realloc(sqia->is, gbuf*sizeof(i_s));
                for(i=gbuf-GBUF;i<gbuf;++i) {
                    sqia->is[i].ibf=GBUF;
                    sqia->is[i].sbf=GBUF;
                    sqia->is[i].id=calloc(sqia->is[i].ibf, sizeof(char));
                    sqia->is[i].sq=calloc(sqia->is[i].sbf, sizeof(char));
                }
            }
            sqia->is[sqidx].idx=sqidx;

            /* resetting stuff */
            sqia->is[sqidx].sylen=0;
            cx=0;
        } else if (IGLINE==1) {
            CONDREALLOC_(cx, sqia->is[sqidx].ibf, GBUF, sqia->is[sqidx].id, char);
            sqia->is[sqidx].id[cx++]=c;
        } else if (IGLINE==0) {
            CONDREALLOC(sqia->is[sqidx].sylen, sqia->is[sqidx].sbf, GBUF, sqia->is[sqidx].sq, char);
            sqia->is[sqidx].sq[sqia->is[sqidx].sylen++]=c;
        }
    }
    fclose(fin);

    /* the last sequence requires special treatment */
    CONDREALLOC_(cx, sqia->is[sqidx].ibf, GBUF, sqia->is[sqidx].id, char);
    sqia->is[sqidx].id[cx]='\0';
    CONDREALLOC_(sqia->is[sqidx].sylen, sqia->is[sqidx].sbf, GBUF, sqia->is[sqidx].sq, char);
    sqia->is[sqidx].sq[sqia->is[sqidx].sylen]='\0';
    sqia->is[sqidx].idz=1+cx;
    sqia->is[sqidx].sqz=1+sqia->is[sqidx].sylen;
    sqia->totbases += sqia->is[sqidx].sylen;
    if(sqia->is[sqidx].sylen > sqia->mx)
        sqia->mx = sqia->is[sqidx].sylen;
    if(sqia->is[sqidx].sylen < sqia->mn)
        sqia->mn = sqia->is[sqidx].sylen;

    /* normalize */
    sqia->numsq=sqidx+1;
    for(i=sqia->numsq;i<gbuf;++i) {
        free(sqia->is[i].id);
        free(sqia->is[i].sq);
    }
    sqia->is=realloc(sqia->is, sqia->numsq*sizeof(i_s));
    return sqia;
}

int main(int argc, char *argv[])
{
    /* argument accounting: remember argc, the number of arguments, _includes_ the executable */
    if(argc!=2) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    int i, j, k;
    i_sa *sqia=faf_to_i_s(argv[1]);

    int uoasz;
    uoa_t *uoa = uniquelens(sqia, &uoasz);

    /* Accumulator strat: run through the uoa, creating a array of arrays of positions */
    unsigned abuf=GBUF;
    int ai=0, acc=0, accsz=0;;
    sia_t *sia=malloc(abuf*sizeof(sia_t));
    for(i=0;i<abuf;++i) {
        sia[i].sisz=0;
        sia[i].tssz=0;
        sia[i].mxssz=0;
        sia[i].sibf=abuf;
        sia[i].si=malloc(sia[i].sibf*sizeof(int));
    }
    unsigned mxsqlen=uoa[uoasz-1].uo; /* the maximum sequence size length */
    unsigned mxasqlen;
    unsigned uoalim=15*uoasz/16;
    printf("%u\n", uoalim);
    for(i=0;i<uoasz;++i) {
        mxasqlen= (i>uoalim)? mxsqlen/10 : 2*mxsqlen ; /* the maximum ALLOWED sequence size length */
        for(j=0;j<uoa[i].uoisz;++j) {
            accsz += uoa[i].uo +MINPADNLEN;
            if(accsz > mxasqlen) {
                accsz -= MINPADNLEN;
                SIACONDREALLOC(ai, abuf, GBUF, sia, k);
                sia[ai].tssz=accsz;
                sia[ai].mxssz = uoa[i].uo;
                accsz=0;
                sia[ai].sisz=acc;
                acc=0;
                ai++;
            }
            CONDREALLOC2(acc, sia[ai].sibf, GBUF, sia[ai].si, int);
            sia[ai].si[acc++] = uoa[i].uoids[j];
        }
    }
    /* last one */
    sia[ai].tssz=mxsqlen;
    sia[ai].mxssz=mxsqlen;
    sia[ai].sisz=1;
    sia[ai].si[0]=uoa[uoasz-1].uoids[uoa[uoasz-1].uoisz-1];

    /* get output file ready: this will hold the stitched output file */
    unsigned fnsz=1+strlen(argv[1]);
    char *tp, *tp2, *foutname=calloc(256, sizeof(char));
    char insertstr[10]="_stitched";
    foutname=realloc(foutname, (fnsz+10)*sizeof(char));
    tp=strrchr(argv[1], '.'); /* the final . position */
    tp2=strrchr(argv[1], '/'); /* the final / position, useful if files are in another directory */
    char fhistname[256]={0}; /* also want to output the histogram */
    if(tp2) {
        sprintf(foutname, "%.*s%s%s", (int)(tp-tp2-1), tp2+1, insertstr, tp);
        sprintf(fhistname, "%.*s%s", (int)(tp-tp2-1), tp2+1, ".uoahist");
    } else { 
        sprintf(foutname, "%.*s%s%s", (int)(tp-argv[1]), argv[1], insertstr, tp);
        sprintf(fhistname, "%.*s%s", (int)(tp-argv[1]), argv[1], ".uoahist");
    }

    printf("%s\n", fhistname); 
    fprtuoahist(fhistname, uoa, uoasz);

    FILE *fout=fopen(foutname, "w");

    int siasz=ai+1;
    for(i=0;i<siasz;++i)
        prtsiaele(sqia, fout, sia[i].si, sia[i].sisz, sia[i].mxssz);

#ifdef DBG
    for(i=0;i<siasz;++i) {
        for(j=0;j<sia[i].sisz;++j) 
            printf("%i ", sia[i].si[j]); 
        printf("... [%u,%u]\n", sia[i].mxssz, sia[i].tssz); 
    }
#endif
    /* free it up free */
    for(i=0;i<abuf;++i)
        free(sia[i].si);
    free(sia);

    fclose(fout);

    /* releasing memory */
    for(j=0; j<uoasz;++j)
        free(uoa[j].uoids);
    free(uoa);

    free_i_sa(&sqia);

    free(foutname);
    return 0;
}
