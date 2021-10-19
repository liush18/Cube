
#include "ppp.h"

#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

static FILE *fpAmb_p=NULL;    /* file pointer of ambiguity */

/* open ambiguity file -------------------------------------------------------*/
extern int logAmbOpen_p(const char *file, const prcopt_t *opt)
{
	if (opt->ionoopt!=IONOOPT_IFLC) return 0;
	createdir(file);
	if (!(fpAmb_p=fopen(file, "w"))) {
		trace(2,"%s: file open error path=%s\n",__func__,file);
		return 0;
	}

	fprintf(fpAmb_p,"%-5s","% SAT");
	fprintf(fpAmb_p,"%14s","WL_AMB(c)");
	fprintf(fpAmb_p,"%14s","WL_N");

	fprintf(fpAmb_p,"%14s","WL_SDAMB(c)");
	fprintf(fpAmb_p,"%14s","STD");
	fprintf(fpAmb_p,"%14s","WL_SDAMB_FIX");

	fprintf(fpAmb_p,"%14s","NL_SDAMB(c)");
	fprintf(fpAmb_p,"%14s","STD");
	fprintf(fpAmb_p,"%14s","NL_SDAMB_FIX");

	fprintf(fpAmb_p,"%14s","IF_AMB(m)");
	fprintf(fpAmb_p,"%14s","STD");
	fprintf(fpAmb_p,"%14s","IF_AMB_FIX");

	fprintf(fpAmb_p,"%14s","AZI(deg)");
	fprintf(fpAmb_p,"%14s","ELE(deg)");
	fprintf(fpAmb_p,"\n");
	fflush(fpAmb_p);

	return 1;
}
/* close ambiguity file ------------------------------------------------------*/
extern void logAmbClose_p(void)
{
	if (fpAmb_p) fclose(fpAmb_p);
	fpAmb_p=NULL;
}
extern void logAmb_p(rtk_t *rtk, const obsd_t *obs, int n)
{
	prcopt_t *opt=&rtk->opt;
	ambc_t *amb;
	ambr_t *ambr;
	ssat_t *ssat;
	double *x=rtk->x,*P=rtk->P,*xa=rtk->xa,*Pa=rtk->Pa,fixed;
	char id[8],time[32];
	int i,k,nx=rtk->nx,sat;

	if (!fpAmb_p) return;

	time2str(obs[0].time,time,2);
	fprintf(fpAmb_p,"%22s%6d%10.2f\n",time,rtk->sol.stat,rtk->sol.ratio);
	for (i=0;i<n&&i<MAXOBS;i++) {
		sat=obs[i].sat;
		satno2id(sat,id);
		k=PIB(sat,1,opt);
		amb=rtk->ambc+sat-1;
		ambr=rtk->ambr+sat-1;
		ssat=rtk->ssat+sat-1;

		fprintf(fpAmb_p,"%5s", id);
		fprintf(fpAmb_p,"%14.4f",amb->LC[0]);
		fprintf(fpAmb_p,"%14d",  amb->n[0]);
		/* wide-lane */
		fprintf(fpAmb_p,"%14.4f",ambr->floamb[0]);
		fprintf(fpAmb_p,"%14.4f",ambr->flostd[0]);
		fprintf(fpAmb_p,"%14d",  ambr->fixamb[0]);
		/* narr-lane */
		fprintf(fpAmb_p,"%14.4f",ambr->floamb[1]);
		fprintf(fpAmb_p,"%14.4f",ambr->flostd[1]);
		fprintf(fpAmb_p,"%14d",  ambr->fixamb[1]);
		/* if ambiguity */
		fprintf(fpAmb_p,"%14.4f",x[k]);
		fprintf(fpAmb_p,"%14.4f",SQRT(P[k+k*nx]));
		if (ambr->fixflag) fixed=xa[k];
		else fixed=0.0;
		fprintf(fpAmb_p,"%14.4f",fixed);

		fprintf(fpAmb_p,"%14.2f",ssat->azel[0]*R2D);
		fprintf(fpAmb_p,"%14.2f",ssat->azel[1]*R2D);
		fprintf(fpAmb_p,"\n");
	}
	fprintf(fpAmb_p,"\n");
	fflush(fpAmb_p);
}