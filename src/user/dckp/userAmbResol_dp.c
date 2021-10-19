/*------------------------------------------------------------------------------
* userAmbResol.c : ppp ambiguity resolution based on dck model
*-----------------------------------------------------------------------------*/

#include "dckp.h"

#define CONST_AMB       0.001       /* constraint to fixed ambiguity */
#define SQR(x)          ((x)*(x))
#define SQRT(x)			((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define ROUND(x)        (int)floor((x)+0.5)

#define MIN_AMB_PAIR	4		/* min ambiguity pairs, 3 not ok!! */

/* fixed solution ------------------------------------------------------------*/
static int fixAmbSol(const rtk_t *rtk, double *x, double *P, const int *sats,
	int n, int f)
{
	const prcopt_t *opt=&rtk->opt;
	double *v,*H,*R;
	const int nx=rtk->nx;
	int i,j,info,stat=1;

	v=zeros(n,1); H=zeros(nx,n); R=zeros(n,n);
	/* integer ambiguity constraint */
	for (i=0;i<n;i++) {
		j=DPIB(sats[i],f+1,opt);
		v[i]=(double)rtk->amb[sats[i]-1].value[f]-x[j];
		H[j+i*nx]=1.0;
		R[i+i*n]=SQR(CONST_AMB);
	}

	/* update states with constraints */
	if ((info=filter(x,P,H,v,R,nx,n))) {
		trace(2,"%s fix %s solution filter error info=%d\n",rtk->cTime,
			f?"wide-lane":"narrow-lane",info);
		stat=0;
	}

	free(v); free(H); free(R);
	return stat;
}
/* fixing ambiguity by lambda ------------------------------------------------
* return : -1 for lambda error, 0 for amb less than 4, 1 for ok -------------*/
static int _ambLambda(rtk_t *rtk, double *x, double *P, const int *sats,
	int n, int f)
{
	double *a,*Qa,*F,test[2];
	const int nx=rtk->nx;
	int i,j,index[MAXOBS]={0},stat=0;

	if (n<MIN_AMB_PAIR) { 
		trace(2,"%s _ambLambda not enough amb pairs n=%d\n",rtk->cTime,n); 
		return 0; 
	}

	/* reset ambiguity tracking status */
	for (i=0;i<MAXSAT;i++) {
		rtk->amb[i].flag[f]=0;
		rtk->amb[i].value[f]=0;
	}

	/* record ambiguity index to be fixed */
	for (i=0;i<n;i++) index[i]=DPIB(sats[i],f+1,&rtk->opt);

	a=mat(n,1);Qa=mat(n,n);F=mat(n,2);
	/* obtain ambiguity values and co-variances */
	for (i=0;i<n;i++) {
		a[i]=x[index[i]];
		for (j=0;j<n;j++) Qa[j+i*n]=P[index[j]+index[i]*nx];
	}

	test[1]=rtk->opt.thresar[1];
	if (!lambdaAR(n,2,a,Qa,F,test)) {
		trace(2,"%s lambda error\n",rtk->cTime);
	}
	else {
		/* varidation by ratio-test */
		if (test[0]<rtk->opt.thresar[0]||test[1]<rtk->opt.thresar[1]) {
			trace(2,"%s varidation error : n=%d ratio=%.2f rate=%.2f%%\n",
				rtk->cTime,n,test[0],test[1]*100.0);
		}
		else {
			stat=1;
			rtk->sol.ratio=(float)test[0];
			rtk->sol.rate =(float)(test[1]*100.0);
			trace(2,"%s varidation ok    : n=%d ratio=%.2f rate=%.2f%%\n",
				rtk->cTime,n,test[0],test[1]*100.0);
			
			for (i=0;i<n;i++) {
				/* AR flag and fixed values */
				rtk->amb[sats[i]-1].flag[f]=1;
				rtk->amb[sats[i]-1].value[f]=ROUND(F[i]);
			}
		}
	}

	free(a);free(Qa);free(F);
	return stat?n:-1;
}
/* fixing narrow-lane ambiguities by lambda ----------------------------------*/
static int fixAmbLambda(rtk_t *rtk, double *x, double *P, int *sats, int n, int f)
{
	double std[MAXOBS]={0};
	int i,j,k,m,stat=0;

	/* record narrow-lane std  */
	for (i=0;i<n;i++) {
		j=DPIB(sats[i],f+1,&rtk->opt);
		std[i]=SQRT(P[j+j*rtk->nx]);
	}

	m=n;
	while (1) {
		j=_ambLambda(rtk,x,P,sats,m,f);
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
		}
		else if (j==0) return 0;
		else break;
	}

	if (j) if (!fixAmbSol(rtk,x,P,sats,m,f)) {
		for (i=0;i<MAXSAT;i++) {
			rtk->amb[i].flag[f]=0;
			rtk->amb[i].value[f]=0;
		}
		return 0;
	}

	return 1;
}
/* fixing wide-lane ambiguities by round -------------------------------------*/
static int fixAmbRound(rtk_t *rtk, double *x, double *P, int *sats, int n, int f)
{
	const prcopt_t *opt=&rtk->opt;
	ambt_t *amb;
	double value;
	const int nx=rtk->nx;
	int i,j,m;

	/* reset wide-lane ambiguity tracking states */
	for (i=0;i<MAXSAT;i++) {
		rtk->amb[i].flag[f]=0;
		rtk->amb[i].value[f]=0;
	}
	/* fixing wide-lane */
	for (i=0;i<n;i++) {
		j=DPIB(sats[i],f+1,&rtk->opt);
		/* not fixed satellites */
		if (SQRT(P[j+j*nx])>opt->thresar[2]||
			fabs(ROUND(x[j])-x[j])>opt->thresar[2]) sats[i]=0;
	}

	/* tracking */
	for (i=0;i<n;i++) {
		if (!sats[i]) continue;
		amb=rtk->amb+sats[i]-1;
		j=DPIB(sats[i],f+1,&rtk->opt);
		value=ROUND(x[j]);

		/* count times fixing to the same value */
		if (amb->track[f]==value) amb->count[f]++;
		else {
			sats[i]=0;
			amb->track[f]=(int)value;
			amb->count[f]=1;
		}
	}
	/* count fixed sats */
	for (i=m=0;i<n;i++) {
		if (!sats[i]) continue;
		sats[m++]=sats[i];
	}
	/* failure */
	if (m<MIN_AMB_PAIR) {
		trace(2,"%s no enough wide-lane n=%d\n",rtk->cTime,m);
		return 0;
	}
	/* update tracking */
	for (i=0;i<m;i++) {
		j=DPIB(sats[i],f+1,&rtk->opt);
		/* ar flag and fixed values */
		rtk->amb[sats[i]-1].flag[f]=1;
		rtk->amb[sats[i]-1].value[f]=ROUND(x[j]);
	}

	if (!fixAmbSol(rtk,x,P,sats,m,f)) {
		for (i=0;i<MAXSAT;i++) {
			rtk->amb[i].flag[f]=0;
			rtk->amb[i].value[f]=0;
		}
		return 0;
	}

	/* fixed wide-lane number */
	return m;
}
/* select ambiguities to fix -------------------------------------------------*/
static int selAmb(const rtk_t *rtk, const obsd_t *obs, int n, const int *exc, 
	const double *azel, const double *P, int *sats)
{
	const prcopt_t *opt=&rtk->opt;
	double elmask;
	const int nx=rtk->nx;
	int i,m;

	elmask=opt->elmaskar>0.0?opt->elmaskar:opt->elmin;
	for (i=m=0;i<n&&i<MAXOBS;i++) {
		if (exc[i]) continue;
		if (azel[i*2+1]<elmask||rtk->datum[obs[i].sat-1]) continue;
		sats[m]=obs[i].sat;
		m++;
	}

	if (m<MIN_AMB_PAIR) {
		trace(2,"%s selAmb not enough ambiguity n=%d\n",rtk->cTime,m);
		return 0;
	}
	return m;
}
/* ambiguity resolution in decoupled clock model -----------------------------*/
extern int userAmbResol_dp(rtk_t *rtk, const obsd_t *obs, int n, const int *exc,
	const double *azel, double *x, double *P)
{
	int m,sats[MAXOBS]={0};

	if (rtk->opt.modear==ARMODE_OFF) return 0;

	/* select ambiguity subset for ar */
	if (!(m=selAmb(rtk,obs,n,exc,azel,P,sats))) return 0;

	if (!(m=fixAmbRound(rtk,x,P,sats,m,1))) return 0;
	if (!fixAmbLambda(rtk,x,P,sats,m,0)) return 0;

	/* 0:no fix, 1:fixed */
	return 1;
}
