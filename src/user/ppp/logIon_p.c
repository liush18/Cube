
#include "ppp.h"

#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

static FILE *fpIon_p=NULL;    /* file pointer of ambiguity */

/* open ionosphere file ------------------------------------------------------*/
extern int logIonOpen_p(const char *file, const prcopt_t *opt)
{
	if (opt->ionoopt==IONOOPT_IFLC) return 0;
	createdir(file);
	if (!(fpIon_p=fopen(file, "w"))) {
		trace(2,"%s: file open error path=%s\n",__func__,file);
		return 0;
	}

	fprintf(fpIon_p,"%-5s","% SAT");
	fprintf(fpIon_p,"%14s","ION(m)");
	fprintf(fpIon_p,"%14s","ION_STD");
	fprintf(fpIon_p,"%14s","ION_FIXED(m)");
	fprintf(fpIon_p,"%14s","ION_STD_FIXED");
	fprintf(fpIon_p,"%14s","AZI(deg)");
	fprintf(fpIon_p,"%14s","ELE(deg)");
	fprintf(fpIon_p,"\n");
	fflush(fpIon_p);

	return 1;
}
/* close ionosphere file -----------------------------------------------------*/
extern void logIonClose_p(void)
{
	if (fpIon_p) fclose(fpIon_p);
	fpIon_p=NULL;
}
extern void logIon_p(rtk_t *rtk, const obsd_t *obs, int n)
{
	prcopt_t *opt=&rtk->opt;
	ssat_t *ssat;
	double *x=rtk->x,*P=rtk->P,*xa=rtk->xa,*Pa=rtk->Pa;
	char id[8];
	int i,k,nx=rtk->nx,sat,stat=rtk->sol.stat;

	if (!fpIon_p) return;

	fprintf(fpIon_p,"%22s%6d%10.2f\n",rtk->cTime,stat,rtk->sol.ratio);
	for (i=0;i<n&&i<MAXOBS;i++) {
		sat=obs[i].sat;
		satno2id(sat,id);
		k=PII(sat,opt);
		ssat=rtk->ssat+sat-1;

		fprintf(fpIon_p,"%5s", id);
		fprintf(fpIon_p,"%14.4f",x[k]);
		fprintf(fpIon_p,"%14.4f",SQRT(P[k+k*nx]));
		fprintf(fpIon_p,"%14.4f",stat==1?xa[k]:0.0);
		fprintf(fpIon_p,"%14.4f",stat==1?Pa[k+k*nx]:0.0);

		fprintf(fpIon_p,"%14.2f",ssat->azel[0]*R2D);
		fprintf(fpIon_p,"%14.2f",ssat->azel[1]*R2D);
		fprintf(fpIon_p,"\n");
	}
	fprintf(fpIon_p,"\n");
	fflush(fpIon_p);
}