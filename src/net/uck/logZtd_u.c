
#include "uck.h"

#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

static FILE *fpZtd=NULL;    /* file pointer of ztd */

/* open ztd file -------------------------------------------------------------*/
extern int openZtdLog_u(const char *file)
{
	createdir(file);
	if (!(fpZtd=fopen(file,"w"))) {
		trace(2,"openZtdLog: file open error path=%s\n",file);
		return 0;
	}
	return 1;
}
/* close ztd file ------------------------------------------------------------*/
extern void closeZtdLog_u(void)
{
	if (fpZtd) fclose(fpZtd);
	fpZtd=NULL;
}
/* output ztd ----------------------------------------------------------------
 * time rcv ztd (gn ge) ztd std (gn_std  ge_std) -----------------------------*/
extern void logZtd_u(net_t *uck, rcv_t *rcv, const int n)
{
	const prcopt_t *opt=&uck->opt;
	double *x,*P;
	int i,j,k,nx=uck->nx;

	if (!fpZtd) return;

	x=uck->x; P=uck->P;
	for (i=0;i<n&&i<MAXRCV;i++) {
		if (rcv[i].outc) continue;
		j=UIT(i+1,opt);
		fprintf(fpZtd,"%-25s%5s",uck->time,rcv[i].sta.name);
		for (k=0;k<UNT(opt);k++) fprintf(fpZtd,"%12.6f",x[j+k]);
		for (k=0;k<UNT(opt);k++) fprintf(fpZtd,"%12.6f",SQRT(P[(j+k)+(j+k)*nx]));
		fprintf(fpZtd,"\n");
	}
	fflush(fpZtd);
}