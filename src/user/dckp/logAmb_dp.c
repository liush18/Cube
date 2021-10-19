
#include "dckp.h"

#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define ROUND(x)        (int)floor((x)+0.5)

static FILE *fpAmb_dp=NULL;    /* file pointer of ambiguity */

/* open ambiguity file -------------------------------------------------------*/
extern int logAmbOpen_dp(const char *file, const prcopt_t *opt)
{
	createdir(file);
	if (!(fpAmb_dp=fopen(file, "w"))) {
		trace(2,"%s: file open error path=%s\n",__func__,file);
		return 0;
	}
	fprintf(fpAmb_dp,"%5s","% SAT");

	fprintf(fpAmb_dp,"%6s","AR");
	fprintf(fpAmb_dp,"%14s","WL_AMB");
	fprintf(fpAmb_dp,"%14s","STD");
	fprintf(fpAmb_dp,"%14s","WL_AMB_FIX");

	fprintf(fpAmb_dp,"%6s","AR");
	fprintf(fpAmb_dp,"%14s","NL_AMB");
	fprintf(fpAmb_dp,"%14s","STD");
	fprintf(fpAmb_dp,"%14s","NL_AMB_FIX");

	fprintf(fpAmb_dp,"%14s", "AZIMUTH");
	fprintf(fpAmb_dp,"%14s", "ELEVATION");
	fprintf(fpAmb_dp,"\n");
	fflush(fpAmb_dp);

	return 1;
}
/* close ambiguity file ------------------------------------------------------*/
extern void logAmbClose_dp(void)
{
	if (fpAmb_dp) fclose(fpAmb_dp);
	fpAmb_dp=NULL;
}
extern void logAmb_dp(rtk_t *rtk, const obsd_t *obs, int n)
{
	prcopt_t *opt=&rtk->opt;
	ambt_t *amb;
	ssat_t *ssat;
	double *x=rtk->x,*P=rtk->P,fixed;
	char id[8];
	int i,nl_k,wl_k,nx=rtk->nx,sat,stat=rtk->sol.stat;

	if (!fpAmb_dp) return;

	fprintf(fpAmb_dp,"%22s%6d%10.2f\n",rtk->cTime,stat,rtk->sol.ratio);
	for (i=0;i<n&&i<MAXOBS;i++) {
		sat=obs[i].sat;
		satno2id(sat,id);
		nl_k=DPIB(sat,1,opt);
		wl_k=DPIB(sat,2,opt);
		amb=rtk->amb+sat-1;
		ssat=rtk->ssat+sat-1;
		/* satllite */
		fprintf(fpAmb_dp,"%5s", id);
		
		/* wide-lane */
		if (amb->flag[1]) fixed=(double)amb->value[1];
		else fixed=0.0;
		fprintf(fpAmb_dp,"%6d",amb->flag[1]);
		fprintf(fpAmb_dp,"%14.4f",x[wl_k]);
		fprintf(fpAmb_dp,"%14.4f",SQRT(P[wl_k+wl_k*nx]));
		fprintf(fpAmb_dp,"%14.4f", fixed);
		/* narr-lane */
		if (amb->flag[0]) fixed=(double)amb->value[0];
		else fixed=0;
		fprintf(fpAmb_dp,"%6d",amb->flag[0]);
		fprintf(fpAmb_dp,"%14.4f",x[nl_k]);
		fprintf(fpAmb_dp,"%14.4f",SQRT(P[nl_k+nl_k*nx]));
		fprintf(fpAmb_dp,"%14.4f",fixed);
		/* azimuth and elevation */
		fprintf(fpAmb_dp,"%14.2f",ssat->azel[0]*R2D);
		fprintf(fpAmb_dp,"%14.2f",ssat->azel[1]*R2D);
		fprintf(fpAmb_dp,"\n");
	}
	fflush(fpAmb_dp);
}