
#include "uck.h"

#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

static FILE *fpUpd=NULL;    /* file pointer of delay */

/* open satellite and receiver delay file ------------------------------------*/
extern int openUPDLog_u(const char *file)
{
	createdir(file);
	if (!(fpUpd=fopen(file, "w"))) {
		trace(2,"openDlyLog: file open error path=%s\n",file);
		return 0;
	}
	return 1;
}
/* close delay file ----------------------------------------------------------*/
extern void closeUPDLog_c(void)
{
	if (fpUpd) fclose(fpUpd);
	fpUpd=NULL;
}
/* output delay --------------------------------------------------------------
* time | rcv/sat | delay | std -----------------------------------------------*/
extern void logUPD_u(net_t *uck, rcv_t *rcv, const int n)
{
	const prcopt_t *opt=&uck->opt;
	double *x,*P;
	char id[8];
	int i,f,k,nx=uck->nx;

	/* no bias for ionosphere free model */
	if (!fpUpd) return;
	if (opt->ionoopt!=IONOOPT_EST) return;

	x=uck->x; P=uck->P;
	/* output receiver delay */
	for (i=0;i<n&&i<MAXRCV;i++) {

		/* skip invalid rcv */
		if (rcv[i].outc) continue;

		fprintf(fpUpd,"%-25s%5s",uck->time,rcv[i].sta.name);
		for (f=0;f<UNF(opt);f++) {
			k=UIRPB(i+1,f+1,opt);
			fprintf(fpUpd,"%12.6f",x[k]);
		}
		for (f=0;f<UNF(opt);f++) {
			k=UIRPB(i+1,f+1,opt);
			fprintf(fpUpd,"%12.6f",SQRT(P[k+k*nx]));
		}
		fprintf(fpUpd,"\n");
	}

	/* output satellite delay */
	for (i=0;i<MAXSAT;i++) {

		/* skip invalid satllite */
		if (uck->outc[i]) continue;

		satno2id(i+1,id);
		fprintf(fpUpd,"%-25s%5s",uck->time,id);
		for (f=0;f<UNF(opt);f++) {
			k=UISPB(i+1,f+1,opt);
			fprintf(fpUpd,"%12.6f",x[k]);
		}
		for (f=0;f<UNF(opt);f++) {
			k=UISPB(i+1,f+1,opt);
			fprintf(fpUpd,"%12.6f",SQRT(P[k+k*nx]));
		}
		fprintf(fpUpd,"\n");
	}
	fflush(fpUpd);
}