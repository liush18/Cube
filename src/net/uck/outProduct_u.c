
#include "uck.h"

#define SQRT(x)    ((x)<0.0||(x)!=(x)?0.0:sqrt(x))

/* sqrt of covariance --------------------------------------------------------*/
static double sqvar(double covar)
{
	return covar<0.0?-sqrt(-covar):sqrt(covar);
}
/* output clock products -----------------------------------------------------*/
extern void outProduct_u(FILE *fp, const net_t *net, const rcv_t *rcv, int n)
{
	const netsol_t *sol=&net->sol;
	double ep[6];
	char id[8],time[32];
	int i,j;

	if (sol->stat==SOLQ_NONE) return;

	time2epoch(sol->time,ep);
	sprintf(time,"%4d %02d %02d %02d %02d%10.6f",
		(int)ep[0],(int)ep[1],(int)ep[2],(int)ep[3],(int)ep[4],ep[5]);

	/* output receiver clock solution */
	for (i=0;i<n&&i<MAXRCV;i++) {
		if (rcv[i].outc) continue;
		j=UIRC(i+1,&net->opt);
		fprintf(fp,"AR %4s %26s%3d   ",rcv[i].sta.name,time,2);
		fprintf(fp,"%19.12E ",net->x[j]*1E-9);
		fprintf(fp,"%19.12E ",SQRT(net->P[j+j*net->nx])*1E-9);
		fprintf(fp,"\n");
	}
	/* output satellite clock solution */
	for (i=0;i<MAXSAT;i++) {
		if (net->outc[i]) continue;
		j=UISC(i+1,&net->opt);
		satno2id(i+1,id);
		fprintf(fp,"AS %-4s %26s%3d   ",id,time,2);
		fprintf(fp,"%19.12E ",net->x[j]*1E-9);
		fprintf(fp,"%19.12E ",SQRT(net->P[j+j*net->nx])*1E-9);
		fprintf(fp,"\n");
	}
	fflush(fp);
}