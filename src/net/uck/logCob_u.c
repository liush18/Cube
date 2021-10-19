
#include "uck.h"

#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

static FILE *fpCob=NULL;    /* file pointer of code bias file */

/* open code bias file -------------------------------------------------------*/
extern int openCobLog_u(const char *file)
{
	createdir(file);
	if (!(fpCob=fopen(file, "w"))) {
		trace(2,"openCobLog: file open error path=%s\n",file);
		return 0;
	}
	return 1;
}
/* close code bias file ------------------------------------------------------*/
extern void closeCobLog_u(void)
{
	if (fpCob) fclose(fpCob);
	fpCob=NULL;
}
/* output code bias of receiver and satellite ----------------------------------
* time rcv/sat cob1 cob2 ... 
* ----------------------------------------------------------------------------*/
extern void logCob_u(net_t *uck, rcv_t *rcv, const int n)
{
	const prcopt_t *opt=&uck->opt;
	double *x,*P;
	char id[8];
	int i,f,k,nx=uck->nx;

	if (!fpCob) return;
	if (!opt->scb&&!opt->rcb) return;

	x=uck->x; P=uck->P;
	if (opt->rcb) for (i=0;i<n&&i<MAXRCV;i++) {
		if (rcv[i].outc) continue;
		fprintf(fpCob,"%-25s%5s",uck->time,rcv[i].sta.name);
		for (f=0;f<UNF(opt);f++) {
			k=UIRCB(i+1,f+1,opt);
			fprintf(fpCob,"%12.4f",x[k]);
		}
		fprintf(fpCob,"\n");
	}
	if (opt->scb) for (i=0;i<MAXSAT;i++) {
		if (uck->outc[i]) continue;
		satno2id(i+1,id);
		fprintf(fpCob,"%-25s%5s",uck->time,id);
		for (f=0;f<UNF(opt);f++) {
			k=UISCB(i+1,f+1,opt);
			fprintf(fpCob,"%12.4f",x[k]);
		}
		fprintf(fpCob,"\n");
	}
	fflush(fpCob);
}