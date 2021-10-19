
#include "dckp.h"

#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

static FILE *fpRes_dp=NULL;    /* file pointer of ambiguity */

/* open residual file --------------------------------------------------------*/
extern int logResOpen_dp(const char *file, const prcopt_t *opt)
{
	createdir(file);
	if (!(fpRes_dp=fopen(file, "w"))) {
		trace(2,"%s: file open error path=%s\n",__func__,file);
		return 0;
	}

	fprintf(fpRes_dp,"%-5s", "% SAT");
	fprintf(fpRes_dp,"%14s%14s%14s","P_IF(m)","L_IF(m)","WL(c)");

	fprintf(fpRes_dp,"%14s","AZI(deg)");
	fprintf(fpRes_dp,"%14s","ELE(deg)");
	fprintf(fpRes_dp,"%14s","OUT");
	fprintf(fpRes_dp,"\n");
	fflush(fpRes_dp);

	return 1;
}
/* close residual file -------------------------------------------------------*/
extern void logResClose_dp(void)
{
	if (fpRes_dp) fclose(fpRes_dp);
	fpRes_dp=NULL;
}
extern void logRes_dp(rtk_t *rtk, const obsd_t *obs, int n)
{
	prcopt_t *opt=&rtk->opt;
	ssat_t *ssat;
	char id[8];
	int i,nx=rtk->nx,sat;

	if (!fpRes_dp) return;

	fprintf(fpRes_dp,"%22s%6d%10.2f\n",rtk->cTime,rtk->sol.stat,rtk->sol.ratio);
	for (i=0;i<n&&i<MAXOBS;i++) {
		sat=obs[i].sat;
		satno2id(sat,id);
		ssat=rtk->ssat+sat-1;

		fprintf(fpRes_dp,"%5s",id);
		fprintf(fpRes_dp,"%14.4f%14.4f%14.4f",ssat->resp[0],ssat->resc[0],
			ssat->resw[0]);

		fprintf(fpRes_dp,"%14.2f",ssat->azel[0]*R2D);
		fprintf(fpRes_dp,"%14.2f",ssat->azel[1]*R2D);
		fprintf(fpRes_dp,"%14d",ssat->outc[0]);
		fprintf(fpRes_dp,"\n");
	}
	fprintf(fpRes_dp,"\n");
	fflush(fpRes_dp);
}