
#include "uck.h"

#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

static FILE *fpAmb=NULL;    /* file pointer of ambiguity */

/* open ambiguity file -------------------------------------------------------*/
extern int openAmbLog_u(const char *file)
{
	createdir(file);
	if (!(fpAmb=fopen(file, "w"))) {
		trace(2,"amb_open: file open error path=%s\n",file);
		return 0;
	}
	return 1;
}
/* close ambiguity file ------------------------------------------------------*/
extern void closeAmbLog_u(void)
{
	if (fpAmb) fclose(fpAmb);
	fpAmb=NULL;
}
/* output float and fixed ambiguity ------------------------------------------*/
extern void logAmb_u(net_t *uck, rcv_t *rcv, const int n)
{
	const prcopt_t *opt=&uck->opt;
	double *x,*P;
	char id[8];
	int i,j,f,k,sat,nx=uck->nx;

	if (!fpAmb) return;

	x=uck->x; P=uck->P;
	for (i=0;i<n&&i<MAXRCV;i++) {
		if (rcv[i].outc) continue;
		for (j=0;j<rcv[i].nobs&&j<MAXOBS;j++) {
			sat=rcv[i].obs[j].sat;
			satno2id(sat,id);
			fprintf(fpAmb,"%-25s%5s%4s",uck->time,rcv[i].sta.name,id);
			for (f=0;f<UNF(opt);f++) {
				k=UIAMB(i+1,sat,f+1,opt);
				fprintf(fpAmb,"%12.4f",x[k]);
			}
			for (f=0;f<UNF(opt);f++) {
				k=UIAMB(i+1,sat,f+1,opt);
				fprintf(fpAmb,"%12.4f",SQRT(P[k+k*nx]));
			}
			fprintf(fpAmb,"%3d\n",rcv[i].datum[sat-1]);
		}
	}
	fflush(fpAmb);
}