/*------------------------------------------------------------------------------
* userAmbResol_uup.c : ppp ambiguity resolution based on uu model
*-----------------------------------------------------------------------------*/

#include "ppp.h"

#define SQR(x)					((x)*(x))
#define SQRT(x)					((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define ROUND(x)				(int)floor((x)+0.5)

#define CONST_AMB				0.001   /* constraint to fixed ambiguity */
#define MIN_SAT_AR				4		/* min ambiguity pairs, 3 not ok!! */
#define MAX_POS_STD_TO_FIX		0.5		/* ar when position std less than 0.5 */
#define MAX_AMB_STD_TO_FIX		1.5		/* ar when amb std less than 1.5 */

/* remove integers and the fractional part is less than 0.5 ------------------*/
static double removeInt(double a)
{
	double b;

	/* remove integers */
	b=a-(int)a;
	if		(b> 0.5) b-=1.0;
	else if (b<-0.5) b+=1.0;

	return b;
}
/* fixed solution ------------------------------------------------------------*/
static int fixAmbSolUU(const rtk_t *rtk, const obsd_t *obs, int n,
	double *x, double *P, const int *amb, int rsat, int irsat, const int *value,
	const double *bias, const double *varBias, int nv)
{
	double *v,*H,*R;
	char str[32];
	const int nx=rtk->nx,nf=rtk->opt.nf;
	int i,f,m,k,rk[NFREQ],info,stat=1;

	if (nv==0) return 0;

	time2str(rtk->sol.time,str,2);
	for (f=0;f<nf;f++) rk[f]=PIB(rsat,f+1,&rtk->opt);

	v=zeros(nv,1); H=zeros(nx,nv); R=zeros(nv,nv);
	/* integer ambiguity constraint */
	for (i=m=0;i<n&&i<MAXOBS;i++) {
		if (!amb[i]||i==irsat) continue;
		
		for (f=0;f<nf;f++) {
			k=PIB(obs[i].sat,f+1,&rtk->opt);
			v[m]=(double)value[f+i*nf]-(x[k]-x[rk[f]]);
			H[k    +m*nx]= 1.0;
			H[rk[f]+m*nx]=-1.0;
			R[m+m*nv]=SQR(CONST_AMB);
			m++;
		}
	}
	if (m!=nv) {
		trace(2,"%s amb dimension error\n",str);
		stat=0;
	}

	/* update states with constraints */
	if (stat&&(info=filter(x,P,H,v,R,nx,nv))) {
		trace(2,"%s fixing solution filter error info=%d\n",str,info);
		stat=0;
	}

	free(v); free(H); free(R);
	return stat;
}
/* AR by lambda --------------------------------------------------------------*/
static int lambdaARUU(rtk_t *rtk, const obsd_t *obs, int n, const double *x, 
	const double *P, int *amb, int rsat, int irsat, int *value, 
	const double *bias, const double *varBias)
{
	const prcopt_t *opt=&rtk->opt;
	double *B1,*N1,*D,*E,*Q,test[2];
	char str[32];
	int i,j,k,m,rk[NFREQ]={0},stat=0;

	time2str(rtk->sol.time,str,2);
	for (i=0;i<rtk->opt.nf;i++) rk[i]=PIB(rsat,i+1,opt);

	k=n*opt->nf;B1=zeros(k,1);D=zeros(rtk->nx,k);
	/* single difference ambiguity */
	for (i=m=0;i<n&&i<MAXOBS;i++) {
		if (!amb[i]||obs[i].sat==rsat) continue;

		for (j=0;j<rtk->opt.nf;j++) {
			k=PIB(obs[i].sat,j+1,opt);
			/* float single difference ambiguity */
			B1[m]=x[k]-x[rk[j]];
			D[ k   +m*rtk->nx]= 1.0;
			D[rk[j]+m*rtk->nx]=-1.0;
			m++;
		}
	}
	if (m/opt->nf<MIN_SAT_AR) {
		free(B1);free(D);
		return 0;
	}

	N1=zeros(m,2); E=mat(m,rtk->nx); Q=mat(m,m);
	/* co-variance of single difference ambiguity */
	matmul("TN",m,rtk->nx,rtk->nx,1.0,D,P,0.0,E);
	matmul("NN",m,m,rtk->nx,1.0,E,D,0.0,Q);

	tracemat(2,Q,m,m,15,8);
	test[1]=opt->thresar[1];
	if (!lambdaAR(m,2,B1,Q,N1,test)) {
		trace(2,"%s lambda error\n",str);
	}
	else {
		tracemat(2,B1,1,m,10,3);
		tracemat(2,N1,1,m,10,3);
		/* varidation by ratio-test */
		if (test[0]<opt->thresar[0]||test[1]<opt->thresar[1]) {
			trace(2,"%s varidation error : n=%d ratio=%.2f rate=%.2f%%\n",
				str,m/opt->nf,test[0],test[1]*100.0);
		}
		else {
			stat=1;
			rtk->sol.ratio=(float)test[0];
			rtk->sol.rate =(float)(test[1]*100.0);
			trace(2,"%s varidation ok    : n=%d ratio=%.2f rate=%.2f%%\n",
				str,m/opt->nf,rtk->sol.ratio,rtk->sol.rate);

			for (i=m=0;i<n&&i<MAXOBS;i++) {
				if (!amb[i]||obs[i].sat==rsat) continue;
				for (j=0;j<opt->nf;j++) {
					value[j+i*opt->nf]=ROUND(N1[m]);
					m++;
				}
			}
		}
	}

	free(B1); free(N1); free(D); free(E); free(Q);
	return stat?m:-1;
}
/* fixing ambiguity ----------------------------------------------------------*/
static int fixAmbUU(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
	const double *azel, double *x, double *P, int *amb)
{
	const prcopt_t *opt=&rtk->opt;
	double std[MAXOBS][NFREQ]={0},stds[MAXOBS]={0},*bias,*varBias;
	int i,j,k,rsat,irsat,value[MAXOBS*NFREQ]={0},stat=0;

	/* record ambiguity std for each frequency */
	bias=mat(rtk->opt.nf,n); varBias=mat(rtk->opt.nf,n);
	for (i=0;i<n&&i<MAXOBS;i++) {
		if (!amb[i]) continue;
		if (!interpol_phasebias(opt,rtk->sol.time,obs[i].sat,nav,
			&bias[i*opt->nf],&varBias[i*opt->nf])) {
			trace(2,"fixAmbUU: no phase bias for sat=%d\n",obs[i].sat);
			amb[i]=0;
			continue;
		}
		for (j=0;j<rtk->opt.nf;j++) {
			k=PIB(obs[i].sat,j+1,&rtk->opt);
			std[i][j]=SQRT(P[k+k*rtk->nx]);
		}
	}

	while (1) {
		/* select satllite with highest elevation */
		for (i=0,irsat=-1;i<n;i++) {
			if (!amb[i]) continue;
			if (irsat<0||azel[i*2+1]>azel[irsat*2+1]) irsat=i;
		}
		rsat=obs[irsat].sat;
		k=lambdaARUU(rtk,obs,n,x,P,amb,rsat,irsat,value,bias,varBias);

		if (k<0) {
			/* del stds of satellite not for AR */
			for (i=0;i<n&&i<MAXOBS;i++) {
				stds[i]=0.0;
				for (j=0;j<rtk->opt.nf;j++) {
					if (!amb[i]) { std[i][j]=0.0; continue; }
					stds[i]+=SQR(std[i][j]);
				}
				stds[i]=SQRT(stds[i]);	/* multi-frequency stds */
			}
			i=d_max(stds,n);
			if (i<0) { k=0; break; }
			amb[i]=0;
		}
		else break;
	}
	if (k) stat=fixAmbSolUU(rtk,obs,n,x,P,amb,rsat,irsat,value,bias,varBias,k);

	free(bias); free(varBias);
	return stat;
}
/* select satellites to AR ---------------------------------------------------*/
static int selSatAR(const rtk_t *rtk, const obsd_t *obs, int n, const int *exc, 
	const double *azel, const double *P, int *amb)
{
	const prcopt_t *opt=&rtk->opt;
	double elmask;
	const int nx=rtk->nx;
	int i,j,k,count,flag;

	elmask=opt->elmaskar>0.0?opt->elmaskar:opt->elmin;
	for (i=count=0;i<n&&i<MAXOBS;i++) {
		if (exc[i]||azel[i*2+1]<elmask) continue;
		for (j=flag=0;j<opt->nf;j++) {
			k=PIB(obs[i].sat,j+1,opt);
			if (SQRT(P[k+k*rtk->nx])>MAX_AMB_STD_TO_FIX) { 
				flag=1;break;
			}
		}
		if (flag) continue;
		amb[i]=1;
		count++;	/* number of satellite to AR */
	}
	if (count<MIN_SAT_AR) {
		trace(2,"selSatAR: %s not enough satellite to AR n=%d\n",
			time_str(rtk->sol.time,2),count);
		return 0;
	}
	return 1;
}
/* ambiguity resolution in ppp -----------------------------------------------*/
extern int userAmbResol_uup(rtk_t *rtk, const obsd_t *obs, int n, 
	const nav_t *nav, const int *exc, const double *azel, double *x, double *P)
{
	double std[3];
	int i,amb[MAXOBS]={0};

	if (rtk->opt.modear==ARMODE_OFF) return 0;

	for (i=0;i<3;i++) std[i]=sqrt(P[i+i*rtk->nx]);
	if (norm(std,3)>MAX_POS_STD_TO_FIX) return 0;

	if (!selSatAR(rtk,obs,n,exc,azel,P,amb)) return 0;
	if (!fixAmbUU(rtk,obs,n,nav,azel,x,P,amb)) return 0;

	/* 0:no fix, 1:fixed */
	return 1;
}
