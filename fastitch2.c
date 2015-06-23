/* fastitch.c stitches up fasta files */
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PDCHAR 'N' /* put this up at top, what the padding character will be */
#define MINPADNLEN 2 /* normally we also want to apply a minimum N-pad so that the short reads can map independently to the stitched contig components */

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

#define CONDREALLOC2(x, b, c, a, t); \
    if((x)==(b)) { \
        (b) += (c); \
        (a)=realloc((a), (b)*sizeof(t)); \
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

typedef struct
{
    i_s *is;
    int numsq;
} i_sa;

void usage(char *executable)
{
    printf("%s, a program to stitch a fragmented fasta file (usually a reference file).\n", executable);
    printf("Usage: %s <one arguments: the (usually) reference fasta filename.\n", executable);
    printf("Explanation:\n");
    printf("\tThis is a second version that groups the sequences which all have the same size.\n");
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

void uprtseq1(i_s *sqi, int sz, FILE *fp, uoa_t *uoa, int *whichuois, int whichuoisz) /* takes the uoa, and prints - one sequence - only those entries corresponding indices in array whichuois */
{
    int i, ii, j, k;
    char paddingchar=PDCHAR;
    int minnlen=MINPADNLEN;
    unsigned mxclen=0;
    for(i=0;i<whichuoisz;++i) /* let's not assume order-by-uo */
        if(uoa[i].uo >mxclen)
            mxclen=uoa[i].uo;
    int npadlen = (mxclen > minnlen)? mxclen : minnlen;

    for(i=0;i<whichuoisz;++i) { /* all these have to go in one sequence */
        for(ii=0;ii<uoa[whichuois[i]].uoisz;++ii) {
            if( (i==0) & (ii==0) ) {
                fprintf(fp, ">");
                for(k=0;k<sqi[uoa[whichuois[i]].uoids[ii]].idz-1;k++)
                    fprintf(fp, "%c", sqi[uoa[whichuois[i]].uoids[ii]].id[k]);
                fprintf(fp, "\n"); 
            }
            for(k=0;k<sqi[uoa[whichuois[i]].uoids[ii]].sqz-1;k++)
                fprintf(fp, "%c", sqi[uoa[whichuois[i]].uoids[ii]].sq[k]);
            if( (i != whichuoisz-1) | (ii != uoa[whichuois[i]].uoisz-1) )
                for(j=0;j<npadlen;++j) 
                    fputc(paddingchar, fp);
        }
    }
    fprintf(fp, "\n"); 
}

void uprtseq2(i_s *sqi, int sz, FILE *fp, uoa_t *uoa, int *whichuois, int whichuoisz) /* takes the uoa, and prints out the individ sequences seperately: in so-called normal manner */
{
    int i, ii, k;

    for(i=0;i<whichuoisz;++i) { /* all these have to go in one sequence */
        for(ii=0;ii<uoa[whichuois[i]].uoisz;++ii) {
            fprintf(fp, ">");
            for(k=0;k<sqi[uoa[whichuois[i]].uoids[ii]].idz-1;k++)
                fprintf(fp, "%c", sqi[uoa[whichuois[i]].uoids[ii]].id[k]);
            fprintf(fp, "\n"); 
            for(k=0;k<sqi[uoa[whichuois[i]].uoids[ii]].sqz-1;k++)
                fprintf(fp, "%c", sqi[uoa[whichuois[i]].uoids[ii]].sq[k]);
            fprintf(fp, "\n"); 
        }
    }
    fprintf(fp, "\n"); 
}

void uprtseq1f(i_s *sqi, int sz, FILE *fp, uoa_t *uoa, int start, int end) /* uses uprtseq1 to print uoa from start sequence continuously to end (inclusively) */
{
    int i;
    int len=end+1-start;
    int *uois=malloc(len*sizeof(int));
    for(i=0;i<len;++i) 
        uois[i]=start+i;

    uprtseq1(sqi, sz, fp, uoa, uois, len);

    free(uois);
}

void uprtseq2f(i_s *sqi, int sz, FILE *fp, uoa_t *uoa, int start, int end) /* uses uprtseq2 to print uoa from a "start" idx to an "end" idx (inclusively) */
{
    int i;
    int len=end+1-start;
    int *uois=malloc(len*sizeof(int));
    for(i=0;i<len;++i) 
        uois[i]=start+i;

    uprtseq2(sqi, sz, fp, uoa, uois, len);

    free(uois);
}

uoa_t *uniquelens(i_s *sqi, unsigned numsq, int *uoasz_)
{
    unsigned char new;
    unsigned i, j;
    unsigned uoabuf=GBUF;
    int uoasz=0;
    uoa_t *uoa=calloc(uoabuf, sizeof(uoa_t));
    for(i=0; i<numsq;++i) {
        new=1;
        for(j=0; j<uoasz;++j) {
            if(uoa[j].uo == sqi[i].sylen) {
                uoa[j].uoisz++;
                uoa[j].uoids=realloc(uoa[j].uoids, uoa[j].uoisz*sizeof(unsigned));
                uoa[j].uoids[uoa[j].uoisz-1] = i;
                // uoa[uoasz-1].uoids[uoa[j].uoisz-1] = sqi[i].idx;
                new=0;
                break;
            }
        }
        if(new) {
            uoasz++;
            CONDREALLOC2(uoasz, uoabuf, GBUF, uoa, uoa_t);
            uoa[uoasz-1].uo = sqi[i].sylen;
            uoa[uoasz-1].uoisz = 1;
            uoa[uoasz-1].uoids=NULL;
            uoa[uoasz-1].uoids=realloc(uoa[j].uoids, uoa[j].uoisz*sizeof(unsigned));
            uoa[uoasz-1].uoids[uoa[j].uoisz-1] = i;
            // uoa[uoasz-1].uoids[uoa[j].uoisz-1] = sqi[i].idx;
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
    int i, j, c, sqidx;
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

    int uoasz;
    uoa_t *uoa = uniquelens(sqi, numsq, &uoasz);

    unsigned mxsqlen=uoa[uoasz-1].uo; /* the maximum sequence size length */
    int cumul=0;
    int pa_s_e[2]={0};
    int pb_s_e[2]={0};
    int pc_s_e[2]={0};
    pc_s_e[1]=uoasz-1;

    


    for(i=0;i<uoasz;++i) {
        if(uoa[i].uo <400)
            prtseq1(sqi, numsq, fout, uoa[i].uo, uoa[i].uoids, uoa[i].uoisz); /* all sequences of this length will be merged into one sequence */
            /* it would actually be nice to spread them over 2 or three actually */
        } else if( (!pa_s_e[1]) & (uoa[i].uo > 0.02*mxsqlen) ) {
            pa_s_e[1]=i;
            if(pa_s_e[1] < uoasz-1)
                pb_s_e[0]=i+1;
            else 
                break;
        }
        if( (pa_s_e[1]) & (cumul >0.2*numsq) ) {
            pb_s_e[1]=i;
            pc_s_e[0]= (pb_s_e[1] < uoasz-1)? i+1 : uoasz -1;
            break;
        }
    }
#ifdef DBG
    printf("PA: %u -> %u ", pa_s_e[0], pa_s_e[1]); 
    printf("PB: %u -> %u ", pb_s_e[0], pb_s_e[1]); 
    printf("PC: %u -> %u\n", pc_s_e[0], pc_s_e[1]); 
#endif

    // prtseq0(sqi, numsq, fout, uoa[1].uoids, uoa[1].uoisz); /* printout to a file named stitched */
    // prtseq1(sqi, numsq, fout, uoa[3].uo, uoa[3].uoids, uoa[3].uoisz);
    // uprtseq1f(sqi, numsq, fout, uoa, 0, 1);
    // uprtseq2f(sqi, numsq, fout, uoa, uoasz-2, uoasz-1);
    fclose(fout);

    /* releasing memory */
    for(j=0; j<uoasz;++j)
        free(uoa[j].uoids);
    free(uoa);

    for(i=0; i<numsq;++i) {
        free(sqi[i].id);
        free(sqi[i].sq);
    }
    free(sqi);

    free(foutname);
    return 0;
}
