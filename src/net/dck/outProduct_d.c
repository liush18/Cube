
#include "dck.h"


#define SQRT(x)    ((x)<0.0||(x)!=(x)?0.0:sqrt(x))


/* sqrt of covariance --------------------------------------------------------*/
static double sqvar(double covar)
{
	return covar<0.0?-sqrt(-covar):sqrt(covar);
}
/* output decoupled clock model solution -------------------------------------*/
extern void outProduct_d(FILE *fp, const net_t *dck, const rcv_t *rcv, int n)
{
	const prcopt_t *opt=&dck->opt;
	const netsol_t *sol=&dck->sol;
	double ep[6],*x,*P,fact;
	char id[4],time[32];
	const int ref=dck->opt.ircv;
	int i,j,k;

	if (sol->stat==SOLQ_NONE) return;
	x=dck->x; P=dck->P;

	time2epoch(sol->time,ep);
	sprintf(time,"%4d%3d%3d%3d%3d%10.6f",
		(int)ep[0],(int)ep[1],(int)ep[2],(int)ep[3],(int)ep[4],ep[5]);

	/* output rcv clock solution */
	for (i=0;i<n&&i<MAXRCV;i++) {
		if (rcv[i].outc) continue;
		
		fprintf(fp,"AR %4s %26s%3d   ",rcv[i].sta.name,time,CLK_NUM(opt)*2);
		for (j=0;j<CLK_NUM(opt);j++) {
			fact=j==2?1.0:1E-9;
			k=DIRC(i+1,j+1,opt);
			fprintf(fp,"%19.12E ",dck->x[k]*fact);
		}
		for (j=0;j<CLK_NUM(opt);j++) {
			fact=j==2?1.0:1E-9;
			k=DIRC(i+1,j+1,opt);
			fprintf(fp,"%19.12E ",SQRT(dck->P[k+k*dck->nx])*fact);
		}
		fprintf(fp,"\n");
	}
	/* output sat clock solution */
	for (i=0;i<MAXSAT;i++) {
		if (dck->outc[i]) continue;
		satno2id(i+1,id);

		fprintf(fp,"AS %-4s %26s%3d   ",id,time,CLK_NUM(opt)*2);
		for (j=0;j<CLK_NUM(opt);j++) {
			fact=j==2?1.0:1E-9;
			k=DISC(i+1,j+1,opt);
			fprintf(fp,"%19.12E ",dck->x[k]*fact);
		}
		for (j=0;j<CLK_NUM(opt);j++) {
			fact=j==2?1.0:1E-9;
			k=DISC(i+1,j+1,opt);
			fprintf(fp,"%19.12E ",SQRT(dck->P[k+k*dck->nx])*fact);
		}
		fprintf(fp,"\n");
	}
	fflush(fp);
}