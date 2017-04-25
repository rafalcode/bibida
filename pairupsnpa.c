/* pairupsnpa: Pair Up SNP Alignment: a way to print out the pairwise matrix of a SNP FASTA alignment */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined DBG || defined DBG2
#define GBUF 4
#else
#define GBUF 128
#endif
#define SSZ 2 /* CG count, first, AT count second, third are the anomalous characters */
#define HISTBUCKETSZ 10
#define HTCOLWIDTH 120
#define OFFSET 1

#define CONDREALLOC(x, b, c, a, t); \
	if((x)==((b)-1)) { \
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

typedef unsigned char boole;
typedef struct /* onefa */
{
	char *id;
	char *sq;
	unsigned idz, sqz;
} onefa;

typedef struct /* i_s; sequence index and number of symbols */
{
	unsigned idx;
	size_t sylen; /* this is the precise symbol length of the sequence */
	size_t sy[SSZ]; /* used to hold counts of symbols */
	float cgp;
	unsigned ambano[2]; /* number of ambiguous symbols (first), then number of anomalous symbols */
	char *id;
	char *sq;
	unsigned ibf, sbf;
	unsigned idz, sqz;
} i_s; /* sequence index and number of symbols */

typedef struct /* ou uov */
{
	unsigned *ua;
	unsigned ub;
	unsigned *oa;
} uo;

void prtusage(void)
{
	printf("pairupsnpa: Pair-Up SNP Alignment: a way to print out the pairwise matrix of a SNP FASTA alignment.\n".);
	printf("Usage: The first argument must be a single FASTA file with a SNP alignment. If no other arguments are given\n");
	printf("the program will output only the number of taxa, the length of the alignment (which must be uniform) and\n");
	printf("all the taxa names themselves.At the same time a full pairwise matrix ending with .tsv will be generated.\n");
	printf("If more than one arguments is given, the extra arguments are expected to be taxon names which should be\n");
	printf("the only ones included in the output tsv. If the final argument, is not a taxon name but the word \"exclude\"\n");
	printf("then all the taxa except those specified will be included in the output tsv file\n"); 
}

int uniquelens(i_s *sqi, unsigned numsq)
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

	if(ai != 1) {
		printf("Sorry, this program is for alignments, all sequences must be the smae length\n"); 
		exit (EXIT_FAILURE);
	}

#if defined DBG || defined DBG2
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

	return sqi[0].sylen; // the first seuence length will do
}

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

void prtsq(i_s *sqisz, int sz)
{
	int i;
	for(i=0;i<sz;++i) {
		printf("%s\n", sqisz[i].id);
		printf("%s\n", sqisz[i].sq);
	}
	return;
}

void prtpwct(i_s *sqisz, int numsq, int *pwa, int nr, int nc, char *spapad)
{
	int i, j, k, mi, mj;

	printf("%*.*s", 12, 12, " "); 
	for(i=0;i<numsq;++i)
		printf("%*.*s ", 12, 12, sqisz[i].id); 
	printf("\n"); 
	mi=0;
	for(i=0;i<nr;++i) {
		mj=nc-i;
		printf("%d) \"%*.*s\": ", i, 12, 12, sqisz[i].id);
		for(k=0;k<i;++k) 
			printf("%s", spapad); 
		for(j=0;j<mj;++j) {
			printf("%*d ", 12, pwa[mi+j]);
		}
		mi+=numsq-i-1; //multiplier for i
		printf("\n"); 
	}
}

void prtpwct2(i_s *sqisz, int numsq, int *pwa, int nr, int nc, char *spapad, char *htmlfn, int oneln)
{
	int i, j, k, mi, mj;
	FILE *fout=fopen(htmlfn, "w");
	fprintf(fout, "<html>\n");
	fprintf(fout, "\t<head>\n");
	fprintf(fout, "\t\t<title>Pairwise SNP Distance Table</title>");
	fprintf(fout, "\t\t<style>\n");
	fprintf(fout, "\t\t\ttable\n");
	fprintf(fout, "\t\t\t{\n");
	// fprintf(fout, "\t\t\t\ttable-layout: fixed;\n");
	// fprintf(fout, "\t\t\t\twidth: 100px;\n");
	fprintf(fout, "\t\t\t\tborder: 1px solid black;\n");
	fprintf(fout, "\t\t\t\tborder-collapse: collapse;\n");
	fprintf(fout, "\t\t\t}\n");
	fprintf(fout, "\t\t\ttd\n");
	fprintf(fout, "\t\t\t{\n");
	fprintf(fout, "\t\t\t\twidth: %dpx;\n", HTCOLWIDTH);
	fprintf(fout, "\t\t\t\tborder: 1px solid black;\n");
	fprintf(fout, "\t\t\t\ttext-align:right;\n");
	fprintf(fout, "\t\t\t}\n");
	fprintf(fout, "\t\t\tth\n");
	fprintf(fout, "\t\t\t{\n");
	fprintf(fout, "\t\t\t\twidth: %dpx;\n", HTCOLWIDTH);
	fprintf(fout, "\t\t\t\tborder: 1px solid black;\n");
	fprintf(fout, "\t\t\t\ttext-align:right;\n");
	fprintf(fout, "\t\t\t}\n");
	fprintf(fout, "\t\t</style>\n");
	fprintf(fout, "\t</head>\n");
	fprintf(fout, "\t<body>\n");
	fprintf(fout, "\t\t<h3>Pairwise SNP Distance Table</h3>\n");
	fprintf(fout, "\t\t<ul>\n");
	fprintf(fout, "\t\t<li>In the first column of this table we have the list of sequence IDs</li>\n");
	fprintf(fout, "\t\t<li>The names in the first row are those sequence which the names in the first column are measured against</li>\n");
	fprintf(fout, "\t\t<li>A sequence is never measured against itself, only against the other sequences in the FASTA alignment file</li>\n");
	fprintf(fout, "\t\t<li>A total of %d SNP sites were detected in the analysis (the being the uniform length of the sequences in the alignment)</li>", oneln);
	fprintf(fout, "\t\t</ul>\n");
	fprintf(fout, "\t\t<table>\n");

	fprintf(fout, "\t\t\t<tr>");
	fprintf(fout, "<td></td>"); // first col of first row empty
	for(i=1;i<numsq;++i)
		fprintf(fout, "<th>%s</th>", sqisz[i].id); 
	fprintf(fout, "</tr>\n");
	mi=0;
	for(i=0;i<nr;++i) {
		fprintf(fout, "\t\t\t<tr>");
		mj=nc-i;
		fprintf(fout, "<td>%s</td>", sqisz[i].id); 
		for(k=0;k<i;++k) 
			fprintf(fout, "<td></td>");
		for(j=0;j<mj;++j)
			fprintf(fout, "<td>%d</td>", pwa[mi+j]);
		mi+=numsq-i-1; //multiplier for i
		fprintf(fout, "</tr>\n");
	}

	fprintf(fout, "\t\t</table>\n");
	fprintf(fout, "\t</body>\n");
	fprintf(fout, "</html>\n");
	fclose(fout);
}

void prtpwct2tsv(i_s *sqisz, int numsq, int *pwa, int nr, int nc, char *spapad, char *tsvfn, int oneln)
{
	int i, j, k, mi, mj;
	FILE *fout=fopen(tsvfn, "w");
	fprintf(fout, "Pairwise arrangement of a Multiple SNP Alignment.\n");
	fprintf(fout, "In the first column of this table we have the list of sequence IDs.\n");
	fprintf(fout, "The names in the first row are those sequence which the names in the first column are measured against.\n");
	fprintf(fout, "A sequence is never measured against itself, only against the other sequences in the FASTA alignment file.\n");
	fprintf(fout, "A total of %d SNP sites were detected in the analysis (this being the uniform (naturally) length of the sequences in the alignment.\n", oneln);
	fprintf(fout, "\n");

	fprintf(fout, "\t"); // first col of first row empty
	for(i=1;i<numsq;++i)
		fprintf(fout, "%s\t", sqisz[i].id); 
	fprintf(fout, "\n");
	mi=0;
	for(i=0;i<nr;++i) {
		mj=nc-i;
		fprintf(fout, "%s\t", sqisz[i].id); 
		for(k=0;k<i;++k) 
			fprintf(fout, "\t");
		for(j=0;j<mj;++j)
			fprintf(fout, "%d\t", pwa[mi+j]);
		mi+=numsq-i-1; //multiplier for i
		fprintf(fout, "\n");
	}

	fclose(fout);
}

void prtpwct2tsv2(i_s *sqisz, int numsq, int *ma, char *tsvfn) /* this version requires the properly rendered symmetric matrix, generates bassic TSV no comments,
																  as it's expected this will be used with R's cmdscale() */
{
	int i, j;
	FILE *fout=fopen(tsvfn, "w");
	fprintf(fout, "\t"); // first col of first row empty
	for(i=0;i<numsq;++i)
		fprintf(fout, (i==numsq-1)? "%s\n":"%s\t", sqisz[i].id); 
	for(i=0;i<numsq;++i) {
		fprintf(fout, "%s\t", sqisz[i].id); 
		for(j=0;j<numsq;++j)
			fprintf(fout, (j==numsq-1)? "%d\n":"%d\t", ma[numsq*i+j]);
	}
	fclose(fout);
}

void mirut(int *m, int n, boole cpdown) /* Upperleft triangular to lower left trianular entries, or viceversa depending on cpdown */
{
	int i, j, piv;

	for(i=1;i<n;++i) {
		piv=OFFSET+i;
		for(j=0;j<piv-1;++j)
			if(cpdown ==1) // copy upper entries down onto lower
				m[n*i+j]=m[n*j+i];
			else if(cpdown ==0) // copy lower entries up onto upper
				m[n*j+i]=m[n*i+j];
	}
	return;
}

int *rendutma(int *pwa, int nr) /* render as upper triangular matrix ... square matrix returned */
{
	int i, j, mi, mj, m;
	int n=nr+1;
	int nc=nr;
	int *ma=calloc(n*n, sizeof(int));

	printf("The values of the pwa are given:\n"); 
	mi=0;
	for(i=0;i<nr;++i) {
		mj=nc-i; // gradually decreasing extent of the column run
		m=(nr+1)*i;
		for(j=0;j<mj;++j) {
			ma[m+i+j+1]=pwa[mi+j];
			printf("%i ", pwa[mi+j]);
		}
		mi+=n-i-1;
	}
	printf("\n"); 
	mirut(ma, n, 1);
	/*
	   printf("Matrix rendition:\n"); 
	   for(i=0;i<n;++i) {
	   for(j=0;j<n;++j)
	   printf("%d ", ma[n*i+j]);
	   printf("\n"); 
	   }
	   */
	return ma;
}

void prtfa2(onefa *fac)
{
	int i;
	printf("SQZ=%d:", fac->sqz);
	for(i=0;i<3;++i) 
		putchar(fac->sq[i]);
	printf("\n"); 
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
		printf("pwamfa: _P_air_W_ise (verb) _A_ _M_ultiple _F_asta SNP _A_lignment.\n");
		printf("Usage: Pls supply 1+ arguments: Each a multi-fasta SNP alignment filenames. \n");
		exit(EXIT_FAILURE);
	}
	/* general declarations */
	FILE *fin;
	char *tp, *htmlfn=calloc(256, sizeof(char));
	char IGLINE, begline;
	size_t lidx, mxsylen, mnsylen;
	unsigned mxamb, mnamb;
	int i, j, k, c, sqidx;
	int gbuf;
	i_s *sqisz=NULL;
	int whatint;
	unsigned numsq, numano;
	onefa fac;
	unsigned ibf=GBUF, sbf=GBUF;
	fac.id=calloc(ibf, sizeof(char));
	fac.sq=calloc(sbf, sizeof(char));
	int ididx=0, ididx0=0;
	int oneln, npwc, *pwa=NULL, nr, nc, mi, mj; // for the PWCT
	char *spapad="    ";
	int *ma;

	if(!(fin=fopen(argv[1], "r")) ) { /*should one check the extension of the fasta file ? */
		printf("Error. Cannot open \"%s\" file.\n", argv[j]);
		exit(EXIT_FAILURE);
	}

	tp=strrchr(argv[1], '.');
	sprintf(htmlfn, "%.*s%s", (int)(tp-argv[1]), argv[1], ".html");
	/* OK I have to clean this up ... obviuous I want a TSV this time but it's too much trouble to set up getopts and friends! */
	sprintf(htmlfn, "%.*s%s", (int)(tp-argv[1]), argv[1], ".tsv");
	IGLINE=0, begline=1;
	lidx=0, mxsylen=0, mnsylen=0XFFFFFFFFFFFFFFFF;
	mxamb=0, mnamb=0xFFFFFFFF;

	sqidx=-1; /* this is slightly dangerous, you need very much to know what you're doing */
	gbuf=GBUF;
	// sqisz=malloc(gbuf*sizeof(i_s));
	sqisz=realloc(sqisz, gbuf*sizeof(i_s));
	for(i=0;i<gbuf;++i) {
		sqisz[i].ibf=GBUF;
		sqisz[i].sbf=GBUF;
		sqisz[i].id=calloc(sqisz[i].ibf, sizeof(char));
		sqisz[i].sq=calloc(sqisz[i].sbf, sizeof(char));
	}
	for(i=gbuf-GBUF;i<gbuf;++i) {
		sqisz[i].ambano[0]=0;
		sqisz[i].ambano[1]=0;
	}
	whatint=0; /* needs explanation */
	ididx=0;
	ididx0=0;

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
				CONDREALLOC(ididx0, sqisz[sqidx].ibf, GBUF, sqisz[sqidx].id, char);
				sqisz[sqidx].id[ididx0]='\0';
				CONDREALLOC(sqisz[sqidx].sylen, sbf, GBUF, fac.sq, char);
				fac.sq[sqisz[sqidx].sylen]='\0';
				CONDREALLOC(sqisz[sqidx].sylen, sqisz[sqidx].sbf, GBUF, sqisz[sqidx].sq, char);
				sqisz[sqidx].sq[sqisz[sqidx].sylen]='\0';
				fac.idz=1+ididx;
				fac.sqz=1+sqisz[sqidx].sylen;
				sqisz[sqidx].idz=1+ididx0;
				sqisz[sqidx].sqz=1+sqisz[sqidx].sylen;
			}

			sqidx++;
			if(sqidx==gbuf) {
				gbuf+=GBUF;
				sqisz=realloc(sqisz, gbuf*sizeof(i_s));
				for(i=gbuf-GBUF;i<gbuf;++i) {
					sqisz[i].ibf=GBUF;
					sqisz[i].sbf=GBUF;
					sqisz[i].id=calloc(sqisz[i].ibf, sizeof(char));
					sqisz[i].sq=calloc(sqisz[i].sbf, sizeof(char));
				}
			}
			sqisz[sqidx].idx=sqidx;

			/* resetting stuff */
			sqisz[sqidx].sylen=0;
			ididx=0;
			ididx0=0;
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
			CONDREALLOC(ididx0, sqisz[sqidx].ibf, GBUF, sqisz[sqidx].id, char);
			sqisz[sqidx].id[ididx0++]=c;
		} else if (IGLINE==0) {
			CONDREALLOC(sqisz[sqidx].sylen, sbf, GBUF, fac.sq, char);
			fac.sq[sqisz[sqidx].sylen]=c;
			CONDREALLOC(sqisz[sqidx].sylen, sqisz[sqidx].sbf, GBUF, sqisz[sqidx].sq, char);
			sqisz[sqidx].sq[sqisz[sqidx].sylen]=c;
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
	CONDREALLOC(ididx0, sqisz[sqidx].ibf, GBUF, sqisz[sqidx].id, char);
	sqisz[sqidx].id[ididx0]='\0';
	CONDREALLOC(sqisz[sqidx].sylen, sbf, GBUF, fac.sq, char);
	fac.sq[sqisz[sqidx].sylen]='\0';
	CONDREALLOC(sqisz[sqidx].sylen, sqisz[sqidx].sbf, GBUF, sqisz[sqidx].sq, char);
	sqisz[sqidx].sq[sqisz[sqidx].sylen]='\0';
	fac.idz=1+ididx;
	fac.sqz=1+sqisz[sqidx].sylen;
	sqisz[sqidx].idz=1+ididx0;
	sqisz[sqidx].sqz=1+sqisz[sqidx].sylen;

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
	// return 0;
	for(i=numsq;i<gbuf;++i) {
		free(sqisz[i].id);
		free(sqisz[i].sq);
	}
	sqisz=realloc(sqisz, numsq*sizeof(i_s));

	/* check for uniform sequence size, necessary for alignments */
	oneln=uniquelens(sqisz, numsq);
	printf("File %s: numseq=%u, uniform seq size at %d\n", argv[1], numsq, oneln);
#ifdef DBG
	prtsq(sqisz, numsq);
#endif

	/* go ahead with pairwise table */
	npwc=numsq*(numsq-1)/2; // wll known, from the maths.
	nr=numsq-1;
	nc=numsq-1;
	pwa=realloc(pwa, npwc*sizeof(int));
	memset(pwa, 0, npwc*sizeof(int));
	for(k=0;k<oneln;++k) {
		mi=0;
		for(i=0;i<nr;++i) {
			mj=nc-i;
			for(j=0;j<mj;++j)
				if(sqisz[i].sq[k]!=sqisz[i+j+1].sq[k]) {// if the characters are not the same, record a SNP difference
					pwa[mi+j]++;
#ifdef DBG2
					if(i!=0)
						printf("diff@%d ", k); // on which colum
#endif
				}
			mi+=numsq-i-1; //multiplier for i
		}
	}
#ifdef DBG2
	printf("\n"); 
#endif
	/* HTML version */
	// prtpwct2(sqisz, numsq, pwa, nr, nc, spapad, htmlfn, oneln);
	// prtpwct2tsv(sqisz, numsq, pwa, nr, nc, spapad, htmlfn, oneln);
	ma=rendutma(pwa, nr);
	prtpwct2tsv2(sqisz, numsq, ma, htmlfn);
	free(ma);

	for(i=0;i<numsq;++i) {
		free(sqisz[i].id);
		free(sqisz[i].sq);
	}
	free(sqisz);
	free(pwa);
	free(htmlfn);
	return 0;
}
