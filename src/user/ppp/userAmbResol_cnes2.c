/*------------------------------------------------------------------------------
* userAmbResol_cnes.c : ppp ambiguity resolution based cnes model
*-----------------------------------------------------------------------------*/

#include "ppp.h"

#define CONST_AMB       0.001       /* constraint to fixed ambiguity */
#define SQR(x)          ((x)*(x))
#define SQRT(x)			((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define ROUND(x)        (int)floor((x)+0.5)
#define MIN(x,y)        ((x)<(y)?(x):(y))

#define MIN_RATIO				2.0		/* ratio threshold */
#define MIN_SUCRATE				0.99	/* success rate threshold */
#define MIN_AMB_PAIR			4		/* min ambiguity pairs, sat difference */
#define MAX_POS_STD_TO_FIX		0.5		/* ar when position std less than 0.5 */
#define MAX_AMB_STD_TO_FIX		1.5		/* ar when amb std less than 1.5 */
#define MAX_AMB_ROUND_STD		0.15	/* max ambiguity std cycle in round ar */
#define MAX_AMB_ROUND_FRAC		0.50	/* max difference from an integer drop cycle */

#define AMBID_IF		1
#define MIN_ARC_GAP     300.0       /* min arc gap (s) */


/* fixed solution ------------------------------------------------------------*/
static int fix_sol(rtk_t *rtk, const int *sats, int n, double *x, double *P)
{
	double *v,*H,*R;
	const int nx=rtk->nx;
	int i,j,info;

	if (n<=0) return 0;

	v=zeros(n,1); H=zeros(nx,n); R=zeros(n,n);
	/* constraints to fixed ambiguities */
	for (i=0;i<n;i++) {
		j=PIB(sats[i],AMBID_IF,&rtk->opt);
		v[i]=(double)rtk->amb[sats[i]-1].value[0]-x[j];
		H[j+i*nx]=1.0;
		R[i+i*n]=SQR(CONST_AMB);
	}
	/* update states with constraints */
	if ((info=filter(x,P,H,v,R,nx,n))) {
		trace(1,"filter error (info=%d)\n",info);
		free(v); free(H); free(R);
		return 0;
	}

	free(v); free(H); free(R);
	return 1;
}
/* fix narrow-lane ambiguity by ILS ------------------------------------------*/
static int fix_amb_ILS(rtk_t *rtk, const int *sats, int n, double *x, double *P)
{
	double *a,*F,*Qa,test[2];
	const int nx=rtk->nx;
	int i,j,stat=0,index[MAXOBS]={0};

	for (i=0;i<n;i++) {
		index[i]=PIB(sats[i],AMBID_IF,&rtk->opt);
	}
	a=mat(n,1);Qa=mat(n,n);F=mat(n,2);
	/* obtain narr-lane ambiguity values and co-variances */
	for (i=0;i<n;i++) {
		a[i]=x[index[i]];
		for (j=0;j<n;j++) Qa[j+i*n]=P[index[j]+index[i]*nx];
	}

	//tracemat(2,Qa,n,n,10,5);
	test[1]=0;rtk->opt.thresar[1];
	/* integer least square */
	if (!lambdaAR(n,2,a,Qa,F,test)) {
		trace(2,"lambda error\n");
	}
	else {
		/* varidation by ratio-test */
		if (test[0]<rtk->opt.thresar[0]) {
			trace(2,"%s varidation error : n=%d ratio=%.2f rate=%.2f%%\n",
				rtk->cTime,n,test[0],test[1]*100.0);
		}
		else {
			stat=1;
			rtk->sol.ratio=(float)test[0];
			rtk->sol.rate =(float)(test[1]*100.0);
			trace(2,"%s varidation ok    : n=%d ratio=%.2f rate=%.2f%%\n",
				rtk->cTime,n,test[0],test[1]*100.0);
		}
	}

	if (stat) stat=fix_sol(rtk,sats,n,x,P);

	free(a);free(Qa);free(F);
	return stat?n:-1;
}
/* fixing narrow-lane ambiguities by lambda ----------------------------------*/
static int fixNarrAmb(rtk_t *rtk, int *sats, int n, double *x, double *P)
{
	double std[MAXOBS]={0};
	int i,j,k,m;

	/* record narrow-lane ambiguity std  */
	for (i=0;i<n;i++) {
		j=PIB(sats[i],1,&rtk->opt);
		std[i]=SQRT(P[j+j*rtk->nx]);
	}

	m=n;
	while (1) {

		j=fix_amb_ILS(rtk,sats,m,x,P);
		if (j<0) {
			i=d_max(std,m);
			if (i<0) return 0;
			sats[i]=0; std[i]=0.0;

			for (i=k=0;i<m;i++) {
				if (!sats[i]) continue;
				sats[k]=sats[i];
				std[k++]=std[i];
			}
			m=k;
			if (m<MIN_AMB_PAIR) {
				trace(2,"%s not enough narr-lane n=%d\n",rtk->cTime,m);
				return 0;
			}
			continue;
		}
		else if (j==0) return 0;
		else break;
	}

	return 1;
}
/* select ambiguities to fix -------------------------------------------------*/
static int selAmb(const rtk_t *rtk, const obsd_t *obs, int n, const int *exc, 
	const double *azel, int *sats)
{
	const prcopt_t *opt=&rtk->opt;
	double elmask;
	int i,sat,m;

	elmask=opt->elmaskar>0.0?opt->elmaskar:opt->elmin;
	for (i=m=0;i<n&&i<MAXOBS;i++) {
		if (exc[i]) continue;
		if (azel[i*2+1]<elmask||rtk->datum[obs[i].sat-1]) continue;

		sat=obs[i].sat;
		/* not fixing narrow-lane ambiguity if wide-lane not fixed */
		if (!rtk->amb[sat-1].flag[1]) continue;
		sats[m]=obs[i].sat;
		m++;
	}

	if (m<MIN_AMB_PAIR) {
		trace(2,"%s not enough ambiguity pairs n=%d\n",rtk->cTime,m);
		return 0;
	}

	return m;
}
/* ambiguity resolution in integer clock model -------------------------------*/
extern int userAmbResol_cnes2(rtk_t *rtk, const obsd_t *obs, int n, 
	const int *exc, const double *azel, double *x, double *P)
{
	int m,sats[MAXOBS]={0};

	if (rtk->opt.modear==ARMODE_OFF) return 0;

	/* select ambiguity subset for ar */
	if (!(m=selAmb(rtk,obs,n,exc,azel,sats))) return 0;

	if (!fixNarrAmb(rtk,sats,m,x,P)) return 0;

	/* 0:no fix, 1:fixed */
	return 1;
}
