/* the proof-of-concept "ring" .. it works! */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>

#define LBUF 64
// #define NUMINOUT 2 /* the number of in out characters, this will vary if you decide to move the sliding window by more than 1 socus at a time */

struct dcn
{
    char c;
    struct dcn *n;
};
typedef struct dcn dcn_t;

struct strchainode
{
	dcn_t *ri;
	unsigned ocs; /*occurences of the same string in this ring, usefull when held in the chained hash array */
	struct strchainode *n;
};
typedef struct strchainode snod;

void simprt(char *bnam, unsigned bnlen, unsigned stepsz, unsigned totwns, float *gcfra, unsigned fullsqlen)
{
    unsigned i;
    printf("Seq %s (namesz %u seqsz %u):\n", bnam, bnlen, fullsqlen); 
	for(i=0;i<totwns;i+=stepsz) /* i=+10, each tenth window */
		printf("%i: %.3f\n", i+1, gcfra[i]); 
    return;
}

void wigprt(char *bnam, unsigned bnlen, unsigned stepsz, unsigned totwns, float *gcfra, unsigned fullsqlen)
{
    unsigned i;

    printf("browser position %s:1-%u\n", bnam, fullsqlen);
    printf("browser hide all\n");
    printf("track type=wiggle_0 name=\"fixedStep\" description=\"fixedStep format\" visibility=full autoScale=off viewLimits=0.0:1.0 color=50,150,255 yLineMark=0.5 yLineOnOff=on priority=10\n");
    printf("fixedStep chrom=%s start=1 step=%u span=%u\n", bnam, stepsz, stepsz);

	for(i=0;i<totwns;i+=stepsz) /* i=+10, each tenth window */
		printf("%.3f\n", gcfra[i]); 
    return;
}

void wigfprt(char *ofnam, char *bnam, unsigned bnlen, unsigned stepsz, unsigned totwns, float *gcfra, unsigned fullsqlen)
{
    unsigned i;
    FILE *ofnamptr=fopen(ofnam, "w");

    fprintf(ofnamptr, "browser position %s:1-%u\n", bnam, fullsqlen);
    fprintf(ofnamptr, "browser hide all\n");
    fprintf(ofnamptr, "track type=wiggle_0 name=\"fixedStep\" description=\"fixedStep format\" visibility=full autoScale=off viewLimits=0.0:1.0 color=50,150,255 yLineMark=0.5 yLineOnOff=on priority=10\n");
    fprintf(ofnamptr, "fixedStep chrom=%s start=1 step=%u span=%u\n", bnam, stepsz, stepsz);

	for(i=0;i<totwns;i+=stepsz) /* i=+10, each tenth window */
		fprintf(ofnamptr, "%.3f\n", gcfra[i]); 
    fclose(ofnamptr);
    return;
}

char *basn(char *nam, unsigned *topoi)
{
	char *poi=strchr(nam, '.');
	*topoi=poi-nam; /* Mne: topoi ... to the point, the period, the dot */
	char *ban=calloc((*topoi)+1, sizeof(char));
	strncpy(ban, nam, *topoi);
	return ban;
}

void prtring0(dcn_t *tai, size_t sqidx, int *gccou, float *gcfra, int idx)
{
	dcn_t *st=tai;
	printf("sqi %zu: ", sqidx);
	do {
		putchar(st->c);
		st=st->n;
	} while (st !=tai);
	putchar(' ');
	printf("gccou: %i,%i gcfra: %.3f\n", gccou[0], gccou[1], gcfra[idx]);
	return;
}

void prtring(dcn_t *tai, size_t sqidx, char *iocs, int *gccou, float *gcfra, int idx)
{
	dcn_t *st=tai;
	printf("sqi %zu: ", sqidx);
	do {
		putchar(st->c);
		st=st->n;
	} while (st !=tai);
	putchar(' ');
	printf("iocs: %c,%c gccou: %i,%i gcfra: %.3f\n", iocs[0], iocs[1], gccou[0], gccou[1], gcfra[idx]);
	return;
}

float modgcscore(int *gccou, char *iocs) /* a two char-element array is the arg here. The first char is the new symbol to add, and th esecodn is the old one which must be deleted */
{
	int i;
	for(i=0;i<2;++i) 
		switch(iocs[i]) { /* the first charcter is the new one, the second is the old */
			case 'C': case 'c': case 'G': case 'g':
				gccou[0] = (i==0)? gccou[0]+1 : gccou[0]-1; break;
			case 'A': case 'a': case 'T': case 't':
				gccou[1] = (i==0)? gccou[1]+1 : gccou[1]-1; break;
			case 'S': case 's':
				gccou[0] = (i==0)? gccou[0]+2 : gccou[0]-2; break;
			case 'W': case 'w':
				gccou[1] = (i==0)? gccou[1]+2 : gccou[1]-2; break;
			case 'M': case 'm': case 'R': case 'r': case 'Y': case 'y': case 'K': case 'k':
				gccou[0] = (i==0)? gccou[0]+1 : gccou[0]-1;
				gccou[1] = (i==0)? gccou[1]+1 : gccou[1]-1; break;
			case 'V': case 'v': case 'B': case 'b':
				gccou[0] = (i==0)? gccou[0]+2 : gccou[0]-2;
				gccou[1] = (i==0)? gccou[1]+1 : gccou[1]-1; break;
			case 'H': case 'h': case 'D': case 'd': 
				gccou[0] = (i==0)? gccou[0]+1 : gccou[0]-1;
				gccou[1] = (i==0)? gccou[1]+2 : gccou[1]-2; break;
			case 'N': case 'n':
				gccou[0] = (i==0)? gccou[0]+2 : gccou[0]-2;
				gccou[1] = (i==0)? gccou[1]+2 : gccou[1]-2; break;
		}
	return (float)gccou[0]/(gccou[0] + gccou[1]);
}

float ringscore(dcn_t *tai, int *gccou)
{
	dcn_t *st=tai;
	do {
		switch(st->c) {
			case 'C': case 'c': case 'G': case 'g':
				gccou[0]++; break;
			case 'A': case 'a': case 'T': case 't':
				gccou[1]++; break;
			case 'S': case 's':
				gccou[0]+=2; break;
			case 'W': case 'w':
				gccou[1]+=2; break;
			case 'M': case 'm': case 'R': case 'r': case 'Y': case 'y': case 'K': case 'k':
				gccou[0]++;
				gccou[1]++; break;
			case 'V': case 'v': case 'B': case 'b':
				gccou[0]+=2;
				gccou[1]++; break;
			case 'H': case 'h': case 'D': case 'd': 
				gccou[0]++;
				gccou[1]+=2; break;
			case 'N': case 'n':
				gccou[0]+=2;
				gccou[1]+=2; break;
		}
		st=st->n;
	} while (st !=tai);
	return (float)gccou[0]/(gccou[0] + gccou[1]);
}

void freering(dcn_t *mou)
{
	dcn_t *st=mou->n, *nxt;
	while (st !=mou) {
		nxt=st->n;
		free(st);
		st=nxt;
	}
	free(mou);
}

char cmpdcn(dcn_t *st1, dcn_t *st2) /* coutn mumber of times st2 is in st1 */
{
	char yes=1;
	dcn_t *st1_=st1, *st2_=st2;
	do {
		if(st1_->c != st2_->c) {
			yes=0;
			break;
		}
		st1_=st1_->n;
		st2_=st2_->n;
	} while (st1_ !=st1);
	return yes;
}

dcn_t *creadcn(size_t ssz) /* create empty ring of size ssz */
{
	int i;
	dcn_t *mou /* mouth with a tendency to eat the tail*/, *tai /* tail */, *ttmp /* type temporary */;

	mou=malloc(sizeof(dcn_t));
	mou->c='\0';
	tai=malloc(sizeof(dcn_t));
	tai->c='\0';
	ttmp=tai;
	for(i=1;i<ssz-1;++i) {
		ttmp->n=malloc(sizeof(dcn_t));
		ttmp->n->c='\0';
		ttmp=ttmp->n; /* with ->nmove on */
	}
	ttmp->n=mou;
	mou->n=tai; /* this is the vital connection, the underwater one */
	return mou;
}

void procthefile(char *filename, unsigned winsz)
{
	unsigned bnlen;
	char *bnam=basn(filename, &bnlen);
	FILE *fin=fopen(filename, "r");
	int c;
    unsigned stepsz=10;

	unsigned lbuf=LBUF;
	float *gcfra=malloc(lbuf*sizeof(float));

	dcn_t *mou=creadcn(winsz);
	size_t cou=0 /* the window index. final value +1 is total number of windows in the sequence */, sqidx=0 /* although primarily used for printing the ring, it is also used to check the first full scorecalc */;

	char IGLINE=0, begline=1;
	char iocs[2]={'\0', '\0'};
	int gccou[2]={0, 0};
	while( ( (c = fgetc(fin)) != EOF) ) {
		if(c == '\n') {
			IGLINE=0;
			begline=1;
		} else if( (begline==1) & (c == '>') ) {
			IGLINE =1;
			begline=0; 
		} else if (IGLINE==0) {
			mou->c = c;
			if(mou->n->c) { /* checking to see if ring is full */
				if(sqidx==winsz-1) {
					if(cou == lbuf) {
						lbuf += LBUF;
						gcfra=realloc(gcfra, lbuf*sizeof(float));
					}
					gcfra[cou]= ringscore(mou, gccou);
					//                   prtring0(mou->n, sqidx, gccou, gcfra, cou);
				} else if (sqidx > winsz-1) {
					iocs[0]=mou->c;
					if(cou == lbuf) {
						lbuf += LBUF;
						gcfra=realloc(gcfra, lbuf*sizeof(float));
					}
					gcfra[cou]= modgcscore(gccou, iocs);
					//                    prtring(mou->n, sqidx, iocs, gccou, gcfra, cou);
				}
				cou++;
			}
			sqidx++;
			mou=mou->n;
			iocs[1]=mou->c;
		}
	}
	fclose(fin);
	gcfra=realloc(gcfra, cou*sizeof(float)); /* normalize */

    char *of=calloc(bnlen+1+3+1, sizeof(char));
    strncpy(of, bnam, bnlen);
    strncat(of, ".wig", 4);
    wigfprt(of, bnam, bnlen, stepsz, cou, gcfra, sqidx);

	freering(mou);
	free(gcfra);
	free(bnam);
	free(of);

	return;
}

int main(int argc, char *argv[])
{
    DIR *dirp;
    struct dirent *direntp;

    if(argc != 4) {
        fprintf( stderr, "Usage: %s <dir_name> <fileext> <windowsize>\n", argv[0]);
        exit(1);
    }

    if ((dirp = opendir(argv[1])) == NULL) {
        perror(argv[1]);
        exit(2);
    }
    unsigned winsz=atoi(argv[3]);

    char *ppoi=NULL;
    short extsz=strlen(argv[2]);
    /* d_types:DT_BLK blockdev; DT_CHR chardev; DT_DIR straight dir; DT_FIFO namedpipe; DT_LNK symlink; DT_REG regular file; DT_SOCK UNIX domain socket; DT_UNKNOWN file type unknown */
    while ((direntp = readdir(dirp)) != NULL)
        if(direntp->d_type == DT_REG) {
            if( (ppoi=strrchr(direntp->d_name, '.')) != NULL)
                if(strncmp((ppoi+1), argv[2], extsz)==0)
                    procthefile(direntp->d_name, winsz);
        }

    closedir(dirp);
    return 0;
}
