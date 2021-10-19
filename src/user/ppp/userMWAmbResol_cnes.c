
#include "ppp.h"

#define SQR(x)					((x)*(x))
#define SQRT(x)					((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define ROUND(x)				(int)floor((x)+0.5)
#define MAX_AMB_ROUND_STD		0.15	/* max ambiguity std cycle in round ar */
#define MAX_AMB_ROUND_FRAC		0.15	/* max difference from an integer drop cycle */
#define CONST_AMB				0.001   /* constraint to fixed ambiguity */
#define CONSFIXED				5		/* continuous fixed times of ambiguity hold */
#define MIN_SAT_SUBSET			4		/* minimum number of ambiguity subset */

/* fixed solution ------------------------------------------------------------*/
static int fixWideAmbSol(const rtk_t *rtk, int n, const int *sats, double *x, 
	double *P)
{
	double *v,*H,*R;
	const int nw=rtk->nw;
	int i,j,info,stat=1;

	/* no fixed wide-lane ambiguity */
	if (n==0) return 0;

	v=zeros(n,1); H=zeros(nw,n); R=zeros(n,n);
	/* integer ambiguity constraint */
	for (i=0;i<n;i++) {
		j=sats[i]-1;
		v[i]=ROUND(x[j])-x[j];
		H[j+i*nw]=1.0;
		R[i+i*n]=SQR(CONST_AMB);
	}

	/* update states with constraints */
	if (info=filter(x,P,H,v,R,nw,n)) {
		trace(2,"%s fixWideAmbSol %s solution filter error info=%d\n",
			rtk->cTime,info);
		stat=0;
	}

	free(v); free(H); free(R);
	return stat;
}
/* fix wide-lane ambiguity ---------------------------------------------------*/
extern int userMWAmbResol_cnes(rtk_t *rtk, const obsd_t *obs, int n, 
	const int *exc, double *x, double *P)
{
	ambt_t *amb;
	const int nw=rtk->nw;
	double value;
	int i,j,m,k,sat,sats[MAXOBS]={0};

	/* reset mw ambiguity fixed flag */
	for (i=0;i<MAXSAT;i++) {
		rtk->amb[i].flag[1]=0;
		rtk->amb[i].value[1]=0;
	}
	/* fix wide-lane ambiguity */
	for (i=m=0;i<n&&i<MAXOBS;i++) {
		sat=obs[i].sat;
		j=sat-1;
		if (exc[i]) continue;
		if (SQRT(P[j+j*nw])>MAX_AMB_ROUND_STD||
			fabs(ROUND(x[j])-x[j])>MAX_AMB_ROUND_FRAC) continue;
		sats[m++]=sat;
	}

	/* tracking */
	for (i=0;i<m;i++) {
		sat=sats[i];
		amb=rtk->amb+sat-1;
		value=ROUND(x[sat-1]);

		/* count times fixing to the same value */
		if (amb->track[1]==value) amb->count[1]++;
		else {
			sats[i]=0;
			amb->track[1]=(int)value;
			amb->count[1]=1;
		}
	}

	/* count fixed sats */
	for (i=k=0;i<m;i++) {
		if (!sats[i]) continue;
		sats[k++]=sats[i];
	}

	if (k<MIN_SAT_SUBSET) {
		trace(2,"%s no enough wide-lane n=%d\n",rtk->cTime,k);
		return 0;
	}

	/* fixed ambiguity */
	for (i=0;i<k;i++) {
		sat=sats[i];
		amb=rtk->amb+sat-1;
		if (amb->count[1]>CONSFIXED) {
			amb->flag[1]|=1;
			amb->value[1]=amb->track[1];
		}
	}

	if (!fixWideAmbSol(rtk,k,sats,x,P)) {
		for (i=0;i<MAXSAT;i++) {
			rtk->amb[i].flag[1]=0;
			rtk->amb[i].value[1]=0;
		}
		return 0;
	}

	return 1;
}