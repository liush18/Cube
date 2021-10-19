
#include "ppp.h"

#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

static FILE *fpZtd_p=NULL;    /* file pointer of ztd */

/* open ztd file -------------------------------------------------------------*/
extern int logZtdOpen_p(const char *file, const prcopt_t *opt)
{
	createdir(file);
	if (!(fpZtd_p=fopen(file,"w"))) {
		trace(2,"%s: file open error path=%s\n",__func__,file);
		return 0;
	}
	fprintf(fpZtd_p,"%-22s%6s","% TIME","SOL");
	fprintf(fpZtd_p,"%14s","ZTD(m)");
	if (PNT(opt)>1) fprintf(fpZtd_p,"%14s%14s","GN","GE");
	fprintf(fpZtd_p,"%14s","ZTD_STD");
	if (PNT(opt)>1) fprintf(fpZtd_p,"%9s%9s","GN_STD","GE_STD");

	fprintf(fpZtd_p,"%14s","ZTD_FIX");
	if (PNT(opt)>1) fprintf(fpZtd_p,"%14s%14s","GN_FIX","GE_FIX");
	fprintf(fpZtd_p,"%14s","ZTD_STD");
	if (PNT(opt)>1) fprintf(fpZtd_p,"%9s%9s","GN_STD","GE_STD");
	fprintf(fpZtd_p,"\n");
	fflush(fpZtd_p);

	return 1;
}
/* close ztd file ------------------------------------------------------------*/
extern void logZtdClose_p(void)
{
	if (fpZtd_p) fclose(fpZtd_p);
	fpZtd_p=NULL;
}
/* output ztd ----------------------------------------------------------------*/
extern void logZtd_p(rtk_t *rtk)
{
	const prcopt_t *opt=&rtk->opt;
	double *x=rtk->x,*P=rtk->P,*xa=rtk->xa,*Pa=rtk->Pa;
	int i,j,nx=rtk->nx,stat=rtk->sol.stat;

	if (!fpZtd_p) return;

	i=PIT(opt);
	fprintf(fpZtd_p,"%22s%6d",rtk->cTime,stat);
	for (j=0;j<PNT(opt);j++) fprintf(fpZtd_p,"%14.4f",x[i+j]);
	for (j=0;j<PNT(opt);j++) fprintf(fpZtd_p,"%14.4f",SQRT(P[(i+j)+(i+j)*nx]));
	for (j=0;j<PNT(opt);j++) fprintf(fpZtd_p,"%14.4f",stat==1?xa[i+j]:0.0);
	for (j=0;j<PNT(opt);j++) fprintf(fpZtd_p,"%14.4f",stat==1?SQRT(Pa[(i+j)+(i+j)*nx]):0.0);
	fprintf(fpZtd_p,"\n");
	fflush(fpZtd_p);
}
