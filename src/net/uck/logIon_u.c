
#include "uck.h"

#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

static FILE *fpIon=NULL;    /* file pointer of ion */

/* open ion file -------------------------------------------------------------*/
extern int openIonLog_u(const char *file)
{
	createdir(file);
	if (!(fpIon=fopen(file, "w"))) {
		trace(2,"openIONLog: file open error path=%s\n",file);
		return 0;
	}
	return 1;
}
/* close ion file ------------------------------------------------------------*/
extern void closeIonLog_u(void)
{
	if (fpIon) fclose(fpIon);
	fpIon=NULL;
}
/* output slant ionospheric delay at L1 --------------------------------------*/
extern void logIon_u(net_t *uck, rcv_t *rcv, const int n)
{
	const prcopt_t *opt=&uck->opt;
	double *x,*P;
	char id[8];
	int i,j,k,sat,nx=uck->nx;

	if (!fpIon) return;
	if (opt->ionoopt!=IONOOPT_EST) return;

	x=uck->x; P=uck->P;
	for (i=0;i<n&&i<MAXRCV;i++) {
		if (rcv[i].outc) continue;
		for (j=0;j<rcv[i].nobs&&j<MAXOBS;j++) {
			sat=rcv[i].obs[j].sat;
			k=UII(i+1,sat,opt); satno2id(sat,id);
			fprintf(fpIon,"%-25s%5s%4s%12.4f%12.4f%12.6f%12.6f",
				uck->time,rcv[i].sta.name,id,rcv[i].ssat[sat-1].azel[0]*R2D,
				rcv[i].ssat[sat-1].azel[1]*R2D,x[k],SQRT(P[k+k*nx]));
			fprintf(fpIon,"%3d\n",rcv[i].datum[sat-1]);
		}
	}
	fflush(fpIon);
}