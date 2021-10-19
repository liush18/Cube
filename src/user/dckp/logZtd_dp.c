
#include "dckp.h"

#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

static FILE *fpZtd_dp=NULL;    /* file pointer of ztd */

/* open ztd file -------------------------------------------------------------*/
extern int logZtdOpen_dp(const char *file, const prcopt_t *opt)
{
	createdir(file);
	if (!(fpZtd_dp=fopen(file,"w"))) {
		trace(2,"%s: file open error path=%s\n",__func__,file);
		return 0;
	}
	fprintf(fpZtd_dp,"%-22s%6s","% TIME","SOL");
	fprintf(fpZtd_dp,"%14s","ZTD");
	if (DPNT(opt)>1) fprintf(fpZtd_dp,"%14s%14s","GN","GE");
	fprintf(fpZtd_dp,"%9s", "ZTD_STD");
	if (DPNT(opt)>1) fprintf(fpZtd_dp,"%9s%9s","GN_STD","GE_STD");

	fprintf(fpZtd_dp,"%14s","ZTD_FIX");
	if (DPNT(opt)>1) fprintf(fpZtd_dp,"%14s%14s","GN_FIX","GE_FIX");
	fprintf(fpZtd_dp,"%9s", "ZTD_STD");
	if (DPNT(opt)>1) fprintf(fpZtd_dp,"%9s%9s","GN_STD","GE_STD");
	fprintf(fpZtd_dp,"\n");
	fflush(fpZtd_dp);

	return 1;
}
/* close ztd file ------------------------------------------------------------*/
extern void logZtdClose_dp(void)
{
	if (fpZtd_dp) fclose(fpZtd_dp);
	fpZtd_dp=NULL;
}
/* output ztd ----------------------------------------------------------------*/
extern void logZtd_dp(rtk_t *rtk)
{
	const prcopt_t *opt=&rtk->opt;
	double *x=rtk->x,*P=rtk->P,*xa=rtk->xa,*Pa=rtk->Pa;
	int i,j,nx=rtk->nx,stat=rtk->sol.stat;

	if (!fpZtd_dp) return;

	i=DPIT(opt);
	fprintf(fpZtd_dp,"%22s%6d",rtk->cTime,stat);
	for (j=0;j<DPNT(opt);j++) fprintf(fpZtd_dp,"%14.4f",x[i+j]);
	for (j=0;j<DPNT(opt);j++) fprintf(fpZtd_dp,"%9.4f",SQRT(P[(i+j)+(i+j)*nx]));
	for (j=0;j<DPNT(opt);j++) fprintf(fpZtd_dp,"%14.4f",stat==1?xa[i+j]:0.0);
	for (j=0;j<DPNT(opt);j++) fprintf(fpZtd_dp,"%9.4f",stat==1?SQRT(Pa[(i+j)+(i+j)*nx]):0.0);
	fprintf(fpZtd_dp,"\n");
	fflush(fpZtd_dp);
}
