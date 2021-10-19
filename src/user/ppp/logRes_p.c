
#include "ppp.h"

#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

static FILE *fpRes_p=NULL;    /* file pointer of ambiguity */

/* open residual file --------------------------------------------------------*/
extern int logResOpen_p(const char *file, const prcopt_t *opt)
{
	char code[14],phase[14],out[14];
	int i;

	createdir(file);
	if (!(fpRes_p=fopen(file, "w"))) {
		trace(2,"%s: file open error path=%s\n",__func__,file);
		return 0;
	}

	fprintf(fpRes_p,"%-5s", "% SAT");
	if (PNF(opt)==1) fprintf(fpRes_p,"%14s%14s","P_IF(m)","L_IF(m)");
	else {
		for (i=0;i<PNF(opt);i++) {
			sprintf(code, "P%d(m)",i+1);
			sprintf(phase,"L%d(m)",i+1);
			fprintf(fpRes_p,"%14s%14s",code,phase);
		}
	}

	fprintf(fpRes_p,"%14s","AZI(deg)");
	fprintf(fpRes_p,"%14s","ELE(deg)");
	for (i=0;i<PNF(opt);i++) {
		sprintf(out,"OUT_%d",i+1);
		fprintf(fpRes_p,"%14s",out);
	}
	fprintf(fpRes_p,"\n");
	fflush(fpRes_p);

	return 1;
}
/* close residual file -------------------------------------------------------*/
extern void logResClose_p(void)
{
	if (fpRes_p) fclose(fpRes_p);
	fpRes_p=NULL;
}
extern void logRes_p(rtk_t *rtk, const obsd_t *obs, int n)
{
	prcopt_t *opt=&rtk->opt;
	ssat_t *ssat;
	char id[8];
	int i,j,nx=rtk->nx,sat;

	if (!fpRes_p) return;

	fprintf(fpRes_p,"%22s%6d%10.2f\n",rtk->cTime,rtk->sol.stat,rtk->sol.ratio);
	for (i=0;i<n&&i<MAXOBS;i++) {
		sat=obs[i].sat;
		satno2id(sat,id);
		ssat=rtk->ssat+sat-1;

		fprintf(fpRes_p,"%5s",id);
		for (j=0;j<PNF(opt);j++) {
			fprintf(fpRes_p,"%14.4f%14.4f",ssat->resp[j],ssat->resc[j]);
		}

		fprintf(fpRes_p,"%14.2f",ssat->azel[0]*R2D);
		fprintf(fpRes_p,"%14.2f",ssat->azel[1]*R2D);
		for (j=0;j<PNF(opt);j++) fprintf(fpRes_p,"%14d",ssat->outc[j]);
		fprintf(fpRes_p,"\n");
	}
	fprintf(fpRes_p,"\n");
	fflush(fpRes_p);
}