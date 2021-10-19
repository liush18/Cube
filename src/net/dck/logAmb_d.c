
#include "dck.h"

#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

static FILE *fpAmb_d=NULL;    /* file pointer of ambiguity */

/* open ambiguity file -------------------------------------------------------*/
extern int openAmbLog_d(const char *file)
{
	createdir(file);
	if (!(fpAmb_d=fopen(file, "w"))) {
		trace(2,"amb_open: file open error path=%s\n",file);
		return 0;
	}
	fprintf(fpAmb_d,"%-22s","% TIME(GPST)");
	fprintf(fpAmb_d,"%5s",  "RCV");
	fprintf(fpAmb_d,"%5s",  "SAT");
	fprintf(fpAmb_d,"%14s", "NARR");
	fprintf(fpAmb_d,"%14s", "WIDE");
	fprintf(fpAmb_d,"%9s",  "STD_N");
	fprintf(fpAmb_d,"%9s",  "STD_W");
	fprintf(fpAmb_d,"%6s",  "DTM");
	fprintf(fpAmb_d,"%6s",  "OUTC");
	fprintf(fpAmb_d,"\n");

	return 1;
}
/* close ambiguity file ------------------------------------------------------*/
extern void closeAmbLog_d(void)
{
	if (fpAmb_d) fclose(fpAmb_d);
	fpAmb_d=NULL;
}
extern void logAmb_d(net_t *dck, rcv_t *rcv, const int n)
{
	const prcopt_t *opt=&dck->opt;
	double *x,*P;
	char id[8];
	int i,j,k1,k2,sat,nx=dck->nx;

	x=dck->x; P=dck->P;
	for (i=0;i<n&&i<MAXRCV;i++) {
		if (rcv[i].outc) continue;
		for (j=0;j<rcv[i].nobs&&j<MAXOBS;j++) {
			sat=rcv[i].obs[j].sat;
			satno2id(sat,id);
			k1=DIAMB(i+1,sat,1,opt); k2=DIAMB(i+1,sat,2,opt);
			fprintf(fpAmb_d,"%22s",dck->time);
			fprintf(fpAmb_d,"%5s",rcv[i].sta.name);
			fprintf(fpAmb_d,"%5s",id);
			fprintf(fpAmb_d,"%14.4f%14.4f",x[k1],x[k2]);
			fprintf(fpAmb_d,"%9.4f%9.4f",SQRT(P[k1+k1*nx]),SQRT(P[k2+k2*nx]));
			fprintf(fpAmb_d,"%6d",rcv[i].datum[sat-1]);
			fprintf(fpAmb_d,"%6d\n",rcv[i].ssat[sat-1].outc[0]);
		}
	}
	fflush(fpAmb_d);
}