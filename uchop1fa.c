/* chop1fa ... chops a file with one fasta according to file to positions */
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <dirent.h> 


#ifdef DBG
#define GBUF 6
#define PDCHAR ' ' /* put this up at top, what the padding character will be */
#define MINPADNLEN 3 /* normally we also want to apply a minimum N-pad so that the short reads can map independently to the stitched contig components */
#else
#define GBUF 128
#define PDCHAR 'N' /* put this up at top, what the padding character will be */
#define MINPADNLEN 300 /* normally we also want to apply a minimum N-pad so that the short reads can map independently to the stitched contig components */
#endif

#define WBUF 8
typedef unsigned char boole;
#define PREFSZ 6
#define SQWRAP 60 /* length at which sequence should be wrapped */
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

typedef struct /* siaa_t: holder for arrays of sia_t's */
{
    sia_t *sia;
    int sz;
} siaa_t;

typedef struct /* wseq_t */
{
    size_t *wln;
    size_t wsbuf;
    size_t quan;
    size_t lbuf;
    size_t numl;
    size_t *wpla; /* words per line array number of words per line */
} wseq_t;

wseq_t *create_wseq_t(size_t initsz)
{
    wseq_t *words=malloc(sizeof(wseq_t));
    words->wsbuf = initsz;
    words->quan = initsz;
    words->wln=calloc(words->wsbuf, sizeof(size_t));
    words->lbuf=WBUF;
    words->numl=0;
    words->wpla=calloc(words->lbuf, sizeof(size_t));
    return words;
}

void free_wseq(wseq_t *wa)
{
    free(wa->wln);
    free(wa->wpla);
    free(wa);
}

char *mktmpd(void)
{
    struct timeval tsecs;
    gettimeofday(&tsecs, NULL);
    char *myt=calloc(14, sizeof(char));
    strncpy(myt, "tmpdir_", 7);
    sprintf(myt+7, "%lu", tsecs.tv_usec);

    DIR *d;
    while((d = opendir(myt)) != NULL) { /* see if a directory witht the same name exists ... very low likelihood though */
        gettimeofday(&tsecs, NULL);
        sprintf(myt+7, "%lu", tsecs.tv_usec);
        closedir(d);
    }
    closedir(d);
    mkdir(myt, S_IRWXU);

    return myt;
}

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
    printf("%s, a program which takes a fiel holding only fasta sequence, and chops2 it up according to array of positions.\n", executable);
    printf("Usage: %s <fasta_file_with_1_seq> <text_file_with_start_and_end_positions>\n", executable);
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
            fprintf(fp, "%s", (whichuoisz != 1)? "_STITCHEDSQ\n" : "\n"); 
        }
        for(k=0;k<sqia->is[whichuois[i]].sqz-1;k++)
            fprintf(fp, "%c", sqia->is[whichuois[i]].sq[k]);
        if(i != whichuoisz-1)
            for(k=0;k<npadlen;++k) 
                fputc(paddingchar, fp);
    }
    fprintf(fp, "\n"); 
}

void chopfa(i_sa *sqa, int *mat, int nr)
{
    int i, j, k, fl /* full lines */, rm /* remainder */;
    char *tmpd=mktmpd();
    FILE *fp;
    char fn[128]={0};
    int df, s, e;
    for(i=0;i<nr;++i) {
        sprintf(fn, "%s/%.*s_%03d.fa", tmpd, PREFSZ, sqa->is[0].id, i);
        fp=fopen(fn,"w");
        e = mat[2*i+1]-1; /* note we're expecting the positions to be 1-indexed */
        s = mat[2*i]-1;
        fprintf(fp, ">%s_%d_%d\n", sqa->is[0].id, s+1, e+1);
        df = e - s;
        rm=df%SQWRAP;
        fl=df/SQWRAP;
        if(fl) {
            for(j=0;j<fl;++j) {
                for(k=0;k<SQWRAP;++k) 
                    fputc(sqa->is[0].sq[s+j*SQWRAP+k], fp);
                fputc('\n', fp);
            }
            if(rm) 
                for(k=0;k<rm;++k) 
                    fputc(sqa->is[0].sq[s+j*SQWRAP+k], fp); /* hoping j holds its value! */
            fputc('\n', fp);
        } else if (rm) { /* sequence was actually shorter than SQWRAP */
            for(k=0;k<rm;++k) 
                fputc(sqa->is[0].sq[s+k], fp);
            fputc('\n', fp);
        }
        fclose(fp);
    }
    printf("Please look inside \"%s\" folder for your chopped up fasta sequence rendered into separate files.\n", tmpd);
    free(tmpd);
}

void uchopfa(i_sa *sqa, int *mat, int nr)
{
    int i, j, k, fl /* full lines */, rm /* remainder */;
    char *tmpd=mktmpd();
    FILE *fp;
    char fn[128]={0};
    int df, s, e;
    int prestart=0;
    int postend=sqa->is[0].sylen;
    for(i=0;i<=nr;++i) {
        s = (i==0)? prestart : mat[2*(i-1)+1]-1;
        e = (i==nr)? postend: mat[2*i]-1;
        sprintf(fn, "%s/%.*s_%03d.fa", tmpd, PREFSZ, sqa->is[0].id, i);
        fp=fopen(fn,"w");
        fprintf(fp, ">%s_%d_%d\n", sqa->is[0].id, s+1, e+1);
        df = e - s;
        rm=df%SQWRAP;
        fl=df/SQWRAP;
        if(fl) {
            for(j=0;j<fl;++j) {
                for(k=0;k<SQWRAP;++k) 
                    fputc(sqa->is[0].sq[s+j*SQWRAP+k], fp);
                fputc('\n', fp);
            }
            if(rm) 
                for(k=0;k<rm;++k) 
                    fputc(sqa->is[0].sq[s+j*SQWRAP+k], fp); /* hoping j holds its value! */
            fputc('\n', fp);
        } else if (rm) { /* sequence was actually shorter than SQWRAP */
            for(k=0;k<rm;++k) 
                fputc(sqa->is[0].sq[s+k], fp);
            fputc('\n', fp);
        }
        fclose(fp);
    }
    printf("Please look inside \"%s\" folder for your chopped up fasta sequence rendered into separate files.\n", tmpd);
    free(tmpd);
}

void fprtuoahist(char *fname, uoa_t *uoa, int uoasz)
{
    int j;
    FILE *fhist=fopen(fname, "w");
    for(j=0; j<uoasz;++j)
        fprintf(fhist, "%i %u\n", uoa[j].uo, uoa[j].uoisz);
    fclose(fhist);
}

siaa_t *g_siaa_t(uoa_t *uoa, int uoasz)
{
    /* Accumulator strat: run through the uoa, creating a array of arrays of positions */
    unsigned abuf=GBUF;
    int i, j, k, ai=0, acc=0, accsz=0;;
    siaa_t *s=malloc(sizeof(siaa_t));
    s->sia=malloc(abuf*sizeof(sia_t));
    for(i=0;i<abuf;++i) {
        s->sia[i].sisz=0;
        s->sia[i].tssz=0;
        s->sia[i].mxssz=0;
        s->sia[i].sibf=abuf;
        s->sia[i].si=malloc(s->sia[i].sibf*sizeof(int));
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
                SIACONDREALLOC(ai, abuf, GBUF, s->sia, k);
                s->sia[ai].tssz=accsz;
                s->sia[ai].mxssz = uoa[i].uo;
                accsz=0;
                s->sia[ai].sisz=acc;
                acc=0;
                ai++;
            }
            CONDREALLOC2(acc, s->sia[ai].sibf, GBUF, s->sia[ai].si, int);
            s->sia[ai].si[acc++] = uoa[i].uoids[j];
        }
    }
    /* last one */
    s->sia[ai].tssz=mxsqlen;
    s->sia[ai].mxssz=mxsqlen;
    s->sia[ai].sisz=1;
    s->sia[ai].si[0]=uoa[uoasz-1].uoids[uoa[uoasz-1].uoisz-1];

    s->sz=ai+1;
    /* normalize */
    for(i=s->sz;i<abuf;++i) 
        free(s->sia[i].si);
    s->sia=realloc(s->sia, s->sz*sizeof(sia_t));

    return s;
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

int *processinpf(char *fname, int *m, int *n)
{
    /* In order to make no assumptions, the file is treated as lines containing the same amount of words each,
     * except for lines starting with #, which are ignored (i.e. comments). These words are checked to make sure they contain only floating number-type
     * characters [0123456789+-.] only, one string variable is icontinually written over and copied into a growing floating point array each time */

    /* declarations */
    FILE *fp=fopen(fname,"r");
    int i;
    size_t couc /*count chars */, couw=0 /* count words */, oldcouw = 0;
    char c;
    boole inword=0;
    wseq_t *wa=create_wseq_t(GBUF);
    size_t bwbuf=WBUF;
    char *bufword=calloc(bwbuf, sizeof(char)); /* this is the string we'll keep overwriting. */

    int *mat=malloc(GBUF*sizeof(int));

    while( (c=fgetc(fp)) != EOF) {
        /*  take care of  */
        if( (c== '\n') | (c == ' ') | (c == '\t') | (c=='#')) {
            if( inword==1) { /* we end a word */
                wa->wln[couw]=couc;
                bufword[couc++]='\0';
                bufword = realloc(bufword, couc*sizeof(char)); /* normalize */
                mat[couw]=atoi(bufword);
                couc=0;
                couw++;
            }
            if(c=='#') {
                while( (c=fgetc(fp)) != '\n') ;
                continue;
            } else if(c=='\n') {
                if(wa->numl == wa->lbuf-1) {
                    wa->lbuf += WBUF;
                    wa->wpla=realloc(wa->wpla, wa->lbuf*sizeof(size_t));
                    memset(wa->wpla+(wa->lbuf-WBUF), 0, WBUF*sizeof(size_t));
                }
                wa->wpla[wa->numl] = couw-oldcouw;
                oldcouw=couw;
                wa->numl++;
            }
            inword=0;
        } else if( (inword==0) && ((c == 0x2B) | (c == 0x2D) | ((c >= 0x30) && (c <= 0x39))) ) { /* deal with first character of new word, + and - also allowed */
            if(couw == wa->wsbuf-1) {
                wa->wsbuf += GBUF;
                wa->wln=realloc(wa->wln, wa->wsbuf*sizeof(size_t));
                mat=realloc(mat, wa->wsbuf*sizeof(int));
                for(i=wa->wsbuf-GBUF;i<wa->wsbuf;++i)
                    wa->wln[i]=0;
            }
            couc=0;
            bwbuf=WBUF;
            bufword=realloc(bufword, bwbuf*sizeof(char)); /* don't bother with memset, it's not necessary */
            bufword[couc++]=c; /* no need to check here, it's the first character */
            inword=1;
        } else if( (c >= 0x30) && (c <= 0x39) ) {
            if(couc == bwbuf-1) { /* the -1 so that we can always add and extra (say 0) when we want */
                bwbuf += WBUF;
                bufword = realloc(bufword, bwbuf*sizeof(char));
            }
            bufword[couc++]=c;
        } else {
            printf("Error. Non-integer character detected. This program is only for reading integers\n"); 
            free_wseq(wa);
            exit(EXIT_FAILURE);
        }

    } /* end of big for statement */
    fclose(fp);
    free(bufword);

    /* normalization stage */
    wa->quan=couw;
    wa->wln = realloc(wa->wln, wa->quan*sizeof(size_t)); /* normalize */
    mat = realloc(mat, wa->quan*sizeof(int)); /* normalize */
    wa->wpla= realloc(wa->wpla, wa->numl*sizeof(size_t));

    *m= wa->numl;
    int k=wa->wpla[0];
    for(i=1;i<wa->numl;++i)
        if(k != wa->wpla[i])
            printf("Warning: Numcols is not uniform at %i words per line on all lines. This file has one with %zu.\n", k, wa->wpla[i]); 
    *n= k; 
    free_wseq(wa);

    return mat;
}

int main(int argc, char *argv[])
{
    /* argument accounting: remember argc, the number of arguments, _includes_ the executable */
    if(argc!=3) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    i_sa *sqia=faf_to_i_s(argv[1]);
    if(sqia->numsq != 1) {
        printf("Error: This program only handles single-sequence fasta files. Bailing out.\n");
        exit(EXIT_FAILURE);
    }

    int nr, nc;
    int *mat=processinpf(argv[2], &nr, &nc);
    if(nc != 2) {
        printf("Error: the positions file should only have 2 columns of integers. Bailing out.\n");
        exit(EXIT_FAILURE);
    }

    uchopfa(sqia, mat, nr);

    free_i_sa(&sqia);
    free(mat);

    return 0;
}
