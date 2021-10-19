
#include "dck.h"

#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

static FILE *fpZtd_d=NULL;    /* file pointer of ztd */

/* open ztd file -------------------------------------------------------------*/
extern int openZtdLog_d(const char *file, int opt)
{
	createdir(file);
	if (!(fpZtd_d=fopen(file,"w"))) {
		trace(2,"ztd_open: file open error path=%s\n",file);
		return 0;
	}

	fprintf(fpZtd_d,"%-22s","% TIME(GPST)");
	fprintf(fpZtd_d,"%5s",  "RCV");
	fprintf(fpZtd_d,"%14s", "ZTD");
	if (opt==TROPOPT_ESTG) {
		fprintf(fpZtd_d,"%14s","GRAD_N");
		fprintf(fpZtd_d,"%14s","GRAD_E");
	}
	fprintf(fpZtd_d,"%9s", "STD_ZTD");
	if (opt==TROPOPT_ESTG) {
		fprintf(fpZtd_d,"%9s","STD_GN");
		fprintf(fpZtd_d,"%9s","STD_GE");
	}
	fprintf(fpZtd_d,"%6s","OUTC");
	fprintf(fpZtd_d,"\n");

	return 1;
}
/* close ztd file ------------------------------------------------------------*/
extern void closeZtdLog_d(void)
{
	if (fpZtd_d) fclose(fpZtd_d);
	fpZtd_d=NULL;
}
/* output ztd ----------------------------------------------------------------*/
extern void logZtd_d(net_t *dck, rcv_t *rcv, const int n)
{
	double *x,*P;
	int i,j,k,nx=dck->nx;

	x=dck->x; P=dck->P;

	for (i=0;i<n&&i<MAXRCV;i++) {
		j=DIT(i+1,&dck->opt);
		fprintf(fpZtd_d,"%22s",dck->time);
		fprintf(fpZtd_d,"%5s",rcv[i].sta.name);
		for (k=0;k<DNT(&dck->opt);k++) {
			fprintf(fpZtd_d,"%14.6f",x[j+k]);
		}
		for (k=0;k<DNT(&dck->opt);k++) {
			fprintf(fpZtd_d,"%9.4f",SQRT(P[(j+k)+(j+k)*nx]));
		}
		fprintf(fpZtd_d,"%6d\n",rcv[i].outc);
	}
	fflush(fpZtd_d);
}