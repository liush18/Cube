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
#define MIN_AMB_PAIR			3		/* min ambiguity pairs, sat difference */
#define MAX_POS_STD_TO_FIX		0.5		/* ar when position std less than 0.5 */
#define MAX_AMB_STD_TO_FIX		1.5		/* ar when amb std less than 1.5 */
#define MAX_AMB_ROUND_STD		0.15	/* max ambiguity std cycle in round ar */
#define MAX_AMB_ROUND_FRAC		0.50	/* max difference from an integer drop cycle */

#define AMBID_IF		1
#define MIN_ARC_GAP     300.0       /* min arc gap (s) */

/* wave length of LC (m) -----------------------------------------------------*/
static double lam_LC(int i, int j, int k)
{
	const double f1=FREQ1,f2=FREQ2,f5=FREQ5;

	return CLIGHT/(i*f1+j*f2+k*f5);
}
/* fixed solution ------------------------------------------------------------*/
static int fix_sol(rtk_t *rtk, int sat1, const int *sat2, const double *NC, 
	int n, double *x, double *P)
{
	double *v,*H,*R;
	const int nx=rtk->nx;
	int i,j,k,info;

	if (n<=0) return 0;

	v=zeros(n,1); H=zeros(nx,n); R=zeros(n,n);
	/* constraints to fixed ambiguities */
	for (i=0;i<n;i++) {
		j=PIB(sat1   ,AMBID_IF,&rtk->opt);
		k=PIB(sat2[i],AMBID_IF,&rtk->opt);
		v[i]=NC[i]-(x[j]-x[k]);
		H[j+i*nx]= 1.0;
		H[k+i*nx]=-1.0;
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
static int fix_amb_ILS(rtk_t *rtk, int sat1, const int *sat2, int n, double *x, 
	double *P, double *std)
{
	double C1,C2,*B1,*N1,*NC,*D,*E,*Q,s[2],lam_NL=lam_LC(1,1,0),lam1,lam2;
	const int nx=rtk->nx;
	int i,j,k,NW,info,stat;

	lam1=CLIGHT/FREQ1; lam2=CLIGHT/FREQ2;

	C1= SQR(lam2)/(SQR(lam2)-SQR(lam1));/* alpha */
	C2=-SQR(lam1)/(SQR(lam2)-SQR(lam1));/* beta */

	B1=zeros(n,1); N1=zeros(n,2); D=zeros(nx,n); E=mat(n,nx);
	Q=mat(n,n); NC=mat(n,1);
	for (i=0;i<n;i++) {

		j=PIB(sat1,   AMBID_IF,&rtk->opt);
		k=PIB(sat2[i],AMBID_IF,&rtk->opt);

		/* float narrow-lane ambiguity (cycle) */
		NW=rtk->ambr[sat2[i]-1].fixamb[0];
		B1[i]=(x[j]-x[k]+C2*lam2*NW)/lam_NL;
		N1[i]=ROUND(B1[i]);
		rtk->ambr[sat2[i]-1].floamb[1]=B1[i];

		/* narrow-lane ambiguity transformation matrix */
		D[j+i*nx]= 1.0/lam_NL;
		D[k+i*nx]=-1.0/lam_NL;
	}
	if (n<MIN_AMB_PAIR) {
		trace(2,"%s not enough narr-lane ambiguity n=%d\n",rtk->cTime,n);
		return 0;
	}

	/* covariance of narrow-lane ambiguities */
	/* Q(m,m)=D'(m,nx)P(nx,nx)D(nx,m) */
	matmul("TN",n,nx,nx,1.0,D,P,0.0,E);
	matmul("NN",n,n, nx,1.0,E,D,0.0,Q);

	for (i=0;i<n;i++) {
		std[i]=SQRT(Q[i+i*n]);
		rtk->ambr[sat2[i]-1].flostd[1]=std[i];
	}

	/* integer least square */
	if ((info=lambda(n,2,B1,Q,N1,s))) {
		trace(2,"lambda error: info=%d\n",info);
		free(B1); free(N1); free(D); free(E); free(Q); free(NC);
		return 0;
	}
	if (s[0]<=0.0) {
		free(B1); free(N1); free(D); free(E); free(Q); free(NC);
		return 0;
	}

	rtk->sol.ratio=(float)(MIN(s[1]/s[0],999.9));

	/* varidation by ratio-test */
	if (rtk->opt.thresar[0]>0.0&&rtk->sol.ratio<rtk->opt.thresar[0]) {
		trace(2,"varidation error: %s n=%2d ratio=%8.3f\n",rtk->cTime,n,rtk->sol.ratio);
		free(B1); free(N1); free(D); free(E); free(Q); free(NC);
		return -1;
	}
	trace(2,"varidation ok: %s n=%2d ratio=%8.3f\n",rtk->cTime,n,rtk->sol.ratio);

	/* narrow-lane to iono-free ambiguity */
	for (i=0;i<n;i++) {
		NW=rtk->ambr[sat2[i]-1].fixamb[0];
		rtk->ambr[sat2[i]-1].fixamb[1]=ROUND(N1[i]);
		rtk->ambr[sat2[i]-1].fixflag=1;
		NC[i]=C1*lam1*N1[i]+C2*lam2*(N1[i]-NW);
	}
	/* fixed solution */
	stat=fix_sol(rtk,sat1,sat2,NC,n,x,P);

	free(B1); free(N1); free(D); free(E); free(Q); free(NC);
	return stat?n:-1;
}
/* fixing narrow-lane ambiguities by lambda ----------------------------------*/
static int fixNarrAmb(rtk_t *rtk, int sat1, int *sat2, int n, double *x, 
	double *P)
{
	double std[MAXOBS]={0};
	int i,j,k;

	while (1) {
		
		j=fix_amb_ILS(rtk,sat1,sat2,n,x,P,std);
		if (j<0) {
			i=d_max(std,n);
			if (i<0) return 0;
			sat2[i]=0; std[i]=0.0;
			traceimat(2,sat2,1,n,3);
			for (i=k=0;i<n;i++) {
				if (!sat2[i]) continue;
				sat2[k]=sat2[i];
				std[k++]=std[i];
			}
			n=k;
			traceimat(2,sat2,1,n,3);
			continue;
		}
		else if (j==0) return 0;
		else break;
	}

	return 1;
}
/* fixing wide-lane ambiguities by round -------------------------------------*/
static int fixWideAmb(rtk_t *rtk, const nav_t *nav, int sat1, int sat2)
{
	ambc_t *amb1,*amb2;
	double vW,lam_WL=lam_LC(1,-1,0),BW;
	int NW;

	amb1=rtk->ambc+sat1-1;
	amb2=rtk->ambc+sat2-1;
	if (!amb1->n[0]||!amb2->n[0]) return 0;

	/* satellite difference wide-lane ambiguity */
	BW=(amb1->LC[0]-amb2->LC[0])/lam_WL+(nav->wlbias[sat1-1]-nav->wlbias[sat2-1]);
	NW=ROUND(BW);

	/* variance of wide-lane ambiguity */
	vW=(amb1->LCv[0]/amb1->n[0]+amb2->LCv[0]/amb2->n[0])/SQR(lam_WL);

	rtk->ambr[sat2-1].floamb[0]=BW;
	rtk->ambr[sat2-1].flostd[0]=SQRT(vW);
	if (fabs(NW-BW)<=rtk->opt.thresar[2]&&SQRT(vW)<=rtk->opt.thresar[2]&&
		conffunc(NW,BW,sqrt(vW))>=rtk->opt.thresar[1]) {
		rtk->ambr[sat2-1].fixamb[0]=NW;
		return 1;
	}
	return 0;
}
/* select ambiguities to fix -------------------------------------------------*/
static int selAmb(const rtk_t *rtk, const obsd_t *obs, int n, const int *exc, 
	const double *azel, const double *P, int *sats, int *sat)
{
	const prcopt_t *opt=&rtk->opt;
	double elmask,el[MAXOBS]={0};
	const int nx=rtk->nx;
	int i,j,count;

	elmask=opt->elmaskar>0.0?opt->elmaskar:opt->elmin;
	for (i=count=0;i<n&&i<MAXOBS;i++) {
		if (exc[i]) continue;
		if (azel[i*2+1]<elmask||rtk->datum[obs[i].sat-1]) continue;

		j=PIB(obs[i].sat,AMBID_IF,opt);

		sats[i]=obs[i].sat;
		el[i]=azel[i*2+1];
		count++;
	}

	if (count<MIN_AMB_PAIR) {
		trace(2,"%s not enough ambiguity pairs n=%d\n",
			time_str(rtk->sol.time,2),count);
		return 0;
	}

	/* reference satellite */
	if ((i=d_max(el,n))<0) return 0;
	*sat=sats[i]; sats[i]=0;

	return 1;
}
/* carrier-phase LC (m) ------------------------------------------------------*/
static double L_LC(int i, int j, int k, const double *L)
{
	const double f1=FREQ1,f2=FREQ2,f5=FREQ5;
	double L1,L2,L5;
	if ((i&&!L[0])||(j&&!L[1])||(k&&!L[2])) {
		return 0.0;
	}
	L1=L[0];
	L2=L[1];
	L5=L[2];
	return (i*f1*L1+j*f2*L2+k*f5*L5)/(i*f1+j*f2+k*f5);
}
/* pseudorange LC (m) --------------------------------------------------------*/
static double P_LC(int i, int j, int k, const double *P)
{
	const double f1=FREQ1,f2=FREQ2,f5=FREQ5;
	double P1,P2,P5;

	if ((i&&!P[0])||(j&&!P[1])||(k&&!P[2])) {
		return 0.0;
	}
	P1=P[0];
	P2=P[1];
	P5=P[2];
	return (i*f1*P1+j*f2*P2+k*f5*P5)/(i*f1+j*f2+k*f5);
}
/* noise variance of LC (m) --------------------------------------------------*/
static double var_LC(int i, int j, int k, double sig)
{
	const double f1=FREQ1,f2=FREQ2,f5=FREQ5;

	return (SQR(i*f1)+SQR(j*f2)+SQR(k*f5))/SQR(i*f1+j*f2+k*f5)*SQR(sig);
}
/* average LC ----------------------------------------------------------------*/
static void average_LC(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
	const double *azel)
{
	ambc_t *amb;
	const prcopt_t *opt=&rtk->opt;
	double LC1,var1,sig;
	double dantr[NFREQ]={0},dants[NFREQ]={0},L[NFREQ],P[NFREQ],Lc,Pc;
	int i,j,sat;

	for (i=0;i<n;i++) {
		sat=obs[i].sat;
		if (azel[1+2*i]<rtk->opt.elmin) continue;

		/* only gps */
		if (satsys(sat,NULL)!=SYS_GPS) continue;

		/* corrected phase and code measurements */
		corr_meas(obs+i,nav,azel+i*2,opt,dantr,dants,0.0,L,P,&Lc,&Pc);

		/* triple-freq carrier and code LC (m) */
		LC1=L_LC(1,-1, 0,L)-P_LC(1,1,0,P);
		sig=sqrt(SQR(rtk->opt.err[1])+SQR(rtk->opt.err[2]/sin(azel[1+2*i])));

		/* measurement noise variance (m) */
		var1=var_LC(1,1,0,sig*rtk->opt.eratio[0]);

		/* reset */
		amb=rtk->ambc+sat-1;
		if (rtk->ssat[sat-1].slip[0]||rtk->ssat[sat-1].slip[1]||
			rtk->ssat[sat-1].slip[2]||amb->n[0]==0.0||
			fabs(timediff(amb->epoch[0],obs[0].time))>MIN_ARC_GAP) {

			amb->n[0]=amb->n[1]=amb->n[2]=0;
			amb->LC[0]=amb->LC[1]=amb->LC[2]=0.0;
			amb->LCv[0]=amb->LCv[1]=amb->LCv[2]=0.0;
			amb->fixcnt=0;
			for (j=0;j<MAXSAT;j++) amb->flags[j]=0;
		}
		/* averaging */
		if (LC1) {
			amb->n[0]+=1;
			amb->LC [0]+=(LC1 -amb->LC [0])/amb->n[0];
			amb->LCv[0]+=(var1-amb->LCv[0])/amb->n[0];
		}
		amb->epoch[0]=obs[0].time;
	}
}
/* ambiguity resolution in integer clock model -------------------------------*/
extern int userAmbResol_cnes(rtk_t *rtk, const obsd_t *obs, int n, 
	const nav_t *nav, const int *exc, const double *azel, double *x, double *P)
{
	double std[3];
	int i,m,sats[MAXOBS]={0},sat;

	if (rtk->opt.modear==ARMODE_OFF) return 0;

	for (i=0;i<MAXSAT;i++) {
		rtk->ambr[i].fixflag=0;
		for (m=0;m<NFREQ;m++) {
			rtk->ambr[i].floamb[m]=0.0;
			rtk->ambr[i].flostd[m]=0.0;
			rtk->ambr[i].fixamb[m]=0;
		}
	}

	average_LC(rtk,obs,n,nav,azel);

	/* select ambiguity subset for ar */
	if (!selAmb(rtk,obs,n,exc,azel,P,sats,&sat)) return 0;

	/* fix wide-lane ambiguity */
	for (i=m=0;i<n&&i<MAXOBS;i++) {
		if (!sats[i]) continue;
		if (fixWideAmb(rtk,nav,sat,sats[i])) {
			sats[m++]=sats[i];
		}
	}

	for (i=0;i<3;i++) std[i]=sqrt(P[i+i*rtk->nx]);
	if (norm(std,3)>MAX_POS_STD_TO_FIX) return 0;
	if (m<MIN_AMB_PAIR) {
		trace(2,"%s not enough wide-lane ambiguity n=%d\n",
			time_str(rtk->sol.time,2),m);
		return 0;
	}

	if (!fixNarrAmb(rtk,sat,sats,m,x,P)) return 0;

	/* 0:no fix, 1:fixed */
	return 1;
}
