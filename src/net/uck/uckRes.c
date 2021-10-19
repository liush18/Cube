
#include "uck.h"

#define SQR(x)				(x*x)
#define CONST_SAT_ZMC		1E-3	/* constraint to ZMC */
#define THRES_REJECT_RATIO	4.0		/* threshold ratio of variance to reject */
#define EFACT_GPS_L5		10.0    /* error factor of GPS/QZS L5 */

/* constraint by predicted values --------------------------------------------*/
static int uckConstraint(const net_t *uck, const rcv_t *rcv, int n, 
	const double *x, const double *P, const int *ix, double *vp, double *HP, 
	double *DP, int flag)
{
	const prcopt_t *opt=&uck->opt;
	const int nx=uck->nx;
	int i,j,k,l,f,sat,*id,nc=uck->nx;
	const int uuar=uck->opt.ionoopt==IONOOPT_EST&&opt->modear>ARMODE_OFF;

	id=imat(nc,1); nc=0;
	openRecord_u(flag);
	record_u(TYPE_TIME,uck->time,"",0,0,0,0,0.0,flag);

	/* satellite code bias constraint if prn>=0.0 else white noise */
	if (opt->scb&&opt->prn[7]>=0.0)
	for (i=0;i<MAXSAT;i++) for (f=0;f<UNF(opt);f++) {
		k=UISCB(i+1,f+1,opt);
		if (ix[k]&&x[k]!=0.0&&P[k*(nx+1)]!=0.0) {
			for (l=0;l<nx;l++) HP[l+nc*nx]=k==l?1.0:0.0;
			vp[nc]=0.0;
			id[nc++]=k; 
		}
	}

	/* record satellite phase bias constraint */
	if (uuar) for (i=0;i<MAXSAT;i++) for (f=0;f<UNF(opt);f++) {
		k=UISPB(i+1,f+1,opt);
		if (ix[k]&&x[k]!=0.0&&P[k+k*nx]!=0.0) {
			record_u(TYPE_SCB,"","",i+1,f+1,0,0,x[k],flag);
			for (l=0;l<nx;l++) HP[l+nc*nx]=k==l?1.0:0.0;
			vp[nc]=0.0;
			id[nc++]=k; 
		}
	}

	/* record constraint receiver by receiver */
	for (i=0;i<n&&i<MAXRCV;i++) {

		if (rcv[i].nobs==0) continue;

		/* receiver code bias constraint if prn>=0.0 else white noise */
		if (opt->rcb&&opt->prn[6]>=0.0)
		for (f=0;f<UNF(opt);f++) {
			k=UIRCB(i+1,f+1,opt);
			if (ix[k]&&x[k]!=0.0&&P[k*(nx+1)]!=0.0) {
				for (l=0;l<nx;l++) HP[l+nc*nx]=k==l?1.0:0.0;
				vp[nc]=0.0;
				id[nc++]=k;
			}
		}

		/* record receiver phase bias constraint */
		if (uuar) for (f=0;f<UNF(opt);f++) {
			k=UIRPB(i+1,f+1,opt);
			if (ix[k]&&x[k]!=0.0&&P[k+k*nx]!=0.0) {
				record_u(TYPE_RCB,"",rcv[i].sta.name,0,f+1,0,0,x[k],flag);
				for (l=0;l<nx;l++) HP[l+nc*nx]=k==l?1.0:0.0;
				vp[nc]=0.0;
				id[nc++]=k; 
			}
		}

		/* record tropospheric delay constraint */
		k=UIT(i+1,opt);
		for (j=0;j<UNT(opt);j++) {
			if (ix[k+j]&&x[k+j]!=0.0&&P[(k+j)*(nx+1)]!=0.0) {
				record_u(TYPE_TRP,"",rcv[i].sta.name,0,0,j,0,x[k+j],flag);
				for (l=0;l<nx;l++) HP[l+nc*nx]=(k+j)==l?1.0:0.0;
				vp[nc]=0.0;
				id[nc++]=k+j; 
			}
		}

		/* record ionosphere and ambiguity constraint info */
		for (j=0;j<rcv[i].nobs&&j<MAXOBS;j++) {

			sat=rcv[i].obs[j].sat;
			if (rcv[i].datum[sat-1]<OBS_USE) continue;

			/* ionosphere, esitmated as white noise if var<0.0 */
			if (opt->ionoopt==IONOOPT_EST&&opt->prn[2]>=0.0) {
				k=UII(i+1,sat,opt); 
				if (ix[k]&&x[k]!=0.0&&P[k+k*nx]!=0.0) {
					record_u(TYPE_ION,"",rcv[i].sta.name,sat,0,0,0,x[k],flag);
					for (l=0;l<nx;l++) HP[l+nc*nx]=k==l?1.0:0.0;
					vp[nc]=0.0;
					id[nc++]=k; 
				}
			}

			/* record ambiguity datum */
			if (rcv[i].datum[sat-1]>OBS_USE) {
				for (f=0;f<UNF(opt);f++) {
					k=UIAMB(i+1,sat,f+1,opt);
					record_u(TYPE_AMB,"",rcv[i].sta.name,sat,f+1,0,
						rcv[i].datum[sat-1],x[k],flag);
				}
				continue;
			}
			/* ambiguity estimated */
			for (f=0;f<UNF(opt);f++) {
				k=UIAMB(i+1,sat,f+1,opt);
				if (ix[k]&&x[k]!=0.0&&P[k+k*nx]!=0.0) {
					record_u(TYPE_AMB,"",rcv[i].sta.name,sat,f+1,0,
						rcv[i].datum[sat-1],x[k],flag);
					for (l=0;l<nx;l++) HP[l+nc*nx]=k==l?1.0:0.0;
					vp[nc]=0.0;
					id[nc++]=k; 
				}
			}
		}
	}

	/* constraint co- and variance */
	for (i=0;i<nc;i++) {
		k=id[i];
		for (j=i;j<nc;j++) {
			f=id[j];
			DP[i+j*nc]=DP[j+i*nc]=P[k+f*nx];
		}
	}

	recordVar_u(DP,nc,flag);
	closeRecord_u(flag);
	free(id);
	return nc;
}
/* datum costraint eliminate rank-deficiency ---------------------------------*/
static int datumConstraint_u(const net_t *uck, const rcv_t *rcv, int n, 
	const int *ix, double *vl, double *HL, double *DL, double *var, int nv, 
	int epFlag)
{
	const prcopt_t *opt=&uck->opt;
	int i,j,f,k,nx=uck->nx,s,equ,count[MAXSAT]={0};
	const int uuar=uck->opt.ionoopt==IONOOPT_EST&&opt->modear>ARMODE_OFF;

	/* zero-mean condition imposed on all satellite clocks */
	if (opt->datumType==DTM_ZMCSAT&&opt->estclk) {
		for (i=0;i<nx;i++) HL[i+nv*nx]=0.0;
		for (i=0;i<MAXSAT;i++) {
			k=UISC(i+1,opt);
			if (ix[k]) HL[k+nv*nx]=1.0;
		}
		vl[nv]=0.0;
		var[nv]=SQR(CONST_SAT_ZMC);
		nv++;
	}

	/* zero-mean condition imposed on all satellite phase bias */
	if (uuar&&opt->datumType==DTM_ZMCSAT) {
		for (f=0;f<UNF(opt);f++) {
			for (i=0;i<nx;i++) HL[i+nv*nx]=0.0;
			for (i=0;i<MAXSAT;i++) {
				k=UISPB(i+1,f+1,opt);
				if (ix[k]) HL[k+nv*nx]=1.0;
			}
			vl[nv]=0.0;
			var[nv]=SQR(CONST_SAT_ZMC);
			nv++;
		}
	}

	if (opt->rcb&&opt->scb&&epFlag) {

		for (i=0;i<MAXSAT*UNF(opt);i++) {
			for (j=0;j<nx;j++) HL[j+(nv+i)*nx]=0.0;
			vl[nv+i]=0.0;
			var[nv+i]=SQR(CONST_SAT_ZMC);
		}

		for (i=0;i<n&&i<MAXRCV;i++) {
			if (rcv[i].outc) continue;
			for (j=0;j<rcv[i].nobs&&j<MAXOBS;j++) {
				s=rcv[i].obs[j].sat;
				if (uck->outc[s-1]) continue;
				for (f=0;f<UNF(opt);f++) {
					k=UIRCB(i+1,f+1,opt);
					equ=(s-1)*UNF(opt)+f;
					HL[k+(nv+equ)*nx]=1.0;
					count[s-1]|=1;
				}
			}
		}

		for (i=0;i<MAXSAT;i++) {
			if (count[i]) for (f=0;f<UNF(opt);f++) nv++;
		}
	}

	/* measurement variance */
	for (i=0;i<nv;i++) for (f=0;f<nv;f++) DL[i+f*nv]=(i==f)?var[i]:0.0;

	return nv;
}
/* measurement error variance ------------------------------------------------*/
static double varerr(int sys, double el, int idx, int type, const prcopt_t *opt)
{
	double fact=1.0,sinel=sin(el);

	if (type==0) fact*=opt->eratio[0];	/* pseudorange */
	fact*=sys==SYS_GLO?EFACT_GLO:(sys==SYS_SBS?EFACT_SBS:EFACT_GPS);

	if (sys==SYS_GPS||sys==SYS_QZS) {
		if (idx==2) fact*=EFACT_GPS_L5; /* GPS/QZS L5 error factor */
	}
	if (opt->ionoopt==IONOOPT_IFLC) fact*=3.0;
	return SQR(fact*opt->err[1])+SQR(fact*opt->err[2]/sinel);
}
/* reject satellite with large and max post-fit residual ---------------------*/
static void postRejSat(const char *strtime, int post, rcv_t *rcv, 
	const int *rcvi, const int *obsi, const int *frqi, double *ve, int ne)
{
	double vmax;
	char id[8];
	int j,sat,maxrcv,maxobs,maxfrq,rej;

	vmax=ve[0]; maxrcv=rcvi[0]; maxobs=obsi[0]; maxfrq=frqi[0]; rej=0;
	for (j=1;j<ne;j++) {
		if (fabs(vmax)>=fabs(ve[j])) continue;
		vmax=ve[j]; maxrcv=rcvi[j]; maxobs=obsi[j]; maxfrq=frqi[j]; rej=j;
	}
	sat=rcv[maxrcv].obs[maxobs].sat; satno2id(sat,id);
	trace(2,"%s outlier (%2d) rejected rcv(%4s) sat(%s) %s%d res=%9.4f el=%4.1f\n",
		strtime,post,rcv[maxrcv].sta.name,id,maxfrq%2?"L":"P",maxfrq/2+1,
		vmax,rcv[maxrcv].ssat[sat-1].azel[1]*R2D);
	rcv[maxrcv].datum[sat-1]=OBS_REJ;
	rcv[maxrcv].ssat[sat-1].rejc[maxfrq/2]++;
	rcv[maxrcv].ssat[sat-1].vsat[maxfrq/2]=0;
	ve[rej]=0;
}
/* get phase bias ------------------------------------------------------------*/
static double getBias(const prcopt_t *opt, const double *x, 
	double *HL, int *ix, int ircv, int sat, int frq, double freq)
{
	const double lam=CLIGHT/freq;
	double rbias,sbias;
	int k;

	if (opt->datumType&&ircv==opt->ircv) rbias=0.0;
	else {
		/* receiver carrier-phase bias */
		k=UIRPB(ircv,frq,opt);
		rbias=lam*x[k];
		HL[k]=lam;
		ix[k]|=1;
	}

	/* satellite carrier-phase bias */
	k=UISPB(sat,frq,opt);
	sbias=-lam*x[k];
	HL[k]=-lam; 
	ix[k]|=1;

	return rbias+sbias;
}
/* get ambiguity -------------------------------------------------------------*/
static double getAmb(const double *x, double *HL, int *ix, int k,
	int datum, double freq, int opt)
{
	if (datum==OBS_USE) ix[k]|=1;
	HL[k]=(opt==IONOOPT_IFLC)?1.0:CLIGHT/freq;
	return HL[k]*x[k];
}
/* get ionosphere delay ------------------------------------------------------*/
static double getIon(const double *x, double *HL, int *ix, int k,
	int carr, double freq)
{
	double C=SQR(FREQ1/freq)*(carr?-1.0:1.0);
	HL[k]=C;
	ix[k]|=1;

	return C*x[k];
}
/* get troposphere delay -----------------------------------------------------*/
static double getTrop(const prcopt_t *opt, double *HL, int *ix, int k, 
	const double *dtdx)
{
	int i;

	for (i=0;i<UNT(opt);i++) {
		HL[k+i]=dtdx[i];
		ix[k+i]|=1;
	}
	return 0.0;
}
/* get receiver code bias drift ----------------------------------------------*/
static double getCodeBias(const prcopt_t *opt, const double *x, double *HL, 
	int *ix, int r, int s, int f, int epFlag, int rout, int sout)
{
	double bias=0.0;
	int k;

	if (opt->rcb) {
		k=UIRCB(r,f,opt);
		if (epFlag&&!rout) { 
			HL[k]=CLIGHT*1E-9;
			ix[k]|=1;
		}
		bias=CLIGHT*1E-9*x[k];
	}
	if (opt->scb) {
		k=UISCB(s,f,opt);
		if (epFlag&&!sout) {
			HL[k]=-CLIGHT*1E-9;
			ix[k]|=1;
		}
		bias-=CLIGHT*1E-9*x[k];
	}

	return bias;
}
/* get clock state values ----------------------------------------------------*/
static double getClock(const prcopt_t *opt, const double *x, double *HL, 
	int *ix, int ircv, int sat, double dts)
{
	double cdtr,cdts;
	int k;

	/* receiver clock (ns) */
	if (opt->datumType==DTM_REFRCV&&ircv==opt->ircv) cdtr=0.0;
	else {
		k=UIRC(ircv,opt);
		cdtr=CLIGHT*1E-9*x[k];
		HL[k]=CLIGHT*1E-9; 
		ix[k]|=1;
	}

	/* satellite clock (ns) */
	if (!opt->estclk) cdts=-CLIGHT*dts;
	else {
		k=UISC(sat,opt);
		cdts=-CLIGHT*1E-9*x[k];
		HL[k]=-CLIGHT*1E-9;
		ix[k]|=1;
	}

	return cdtr+cdts;
}
/* calculate residuals -------------------------------------------------------*/
extern int uckRes(int post, const net_t *uck, rcv_t *rcv, int n, 
	const nav_t *nav, const double *x, const double *P, int *ix, double *vl, 
	double *HL, double *DL, double *vp, double *HP, double *DP, int *nc, 
	int flag)
{
	const prcopt_t *opt=&uck->opt;
	ssat_t *ssat;
	double freq,*var,*ve,trpx[3]={0},dtdx[3],dtrp,vart;
	const int uuar=uck->opt.ionoopt==IONOOPT_EST&&opt->modear>ARMODE_OFF;
	int i,j,k,l,idx,carr,sat,sys,stat=1,nv=0,ne=0,nx=uck->nx,*rcvi,*obsi,*frqi;

	/* number of observation equations */
	k=UNF(opt)*2*uck->nobs+1+UNF(opt)+MAXSAT; 
	rcvi=imat(k,1); obsi=imat(k,1); frqi=imat(k,1); var=mat(k,1); ve=mat(k,1);

	for (i=0;i<nx;i++) ix[i]=0;
	for (i=0;i<n&&i<MAXRCV;i++) for (j=0;j<MAXSAT;j++) {
		for (k=0;k<UNF(opt);k++) rcv[i].ssat[j].vsat[k]=0;
	}

	/* receiver by receiver */
	for (i=0;i<n&&i<MAXRCV;i++) for(j=0;j<rcv[i].nobs&&j<MAXOBS;j++) {

		sat=rcv[i].obs[j].sat;
		ssat=&rcv[i].ssat[sat-1];
		sys=satsys(sat,NULL);

		/* skip rejected satellite */
		if (rcv[i].datum[sat-1]<OBS_USE) continue;
		
		matcpy(trpx,x+UIT(i+1,opt),UNT(opt),1);
		tropModel(rcv[i].obs[j].time,rcv[i].pos,ssat->azel,trpx,dtdx,&dtrp,&vart);

		/* frequency by frequency, P1, L1, P2, L2, ... */
		for (k=0;k<2*UNF(opt);k++) {

			idx=k/2;  /* frequency, from 0 */
			carr=k%2; /* obs type, 0 for pseudorange, 1 for carrier-phase */

			/* frequency */
			if ((freq=sat2freq(sat,rcv[i].obs[j].code[idx],nav))==0.0) continue;
			/* observables */
			if (opt->ionoopt==IONOOPT_IFLC) {
				if ((vl[nv]=carr?rcv[i].Lc[j]:rcv[i].Pc[j])==0.0) continue;
			}
			else {
				if ((vl[nv]=carr?rcv[i].L[j][idx]:rcv[i].P[j][idx])==0.0) 
					continue;
			}
			
			/* initialization design matrix */
			for (l=0;l<nx;l++) HL[l+nv*nx]=0.0;
			/* receiver-satellite distance */
			vl[nv]-=rcv[i].r[j]+dtrp;
			/* receiver and satellite clock */
			vl[nv]-=getClock(opt,x,&HL[nv*nx],ix,i+1,sat,rcv[i].dts[j*2]);
			/* receiver code bias drift */
			if (!carr) {
				vl[nv]-=getCodeBias(opt,x,&HL[nv*nx],ix,i+1,sat,idx+1,flag,
					rcv[i].outc,uck->outc[sat-1]);
			}
			/* tropospheric delay (m) */
			vl[nv]-=getTrop(opt,&HL[nv*nx],ix,UIT(i+1,opt),dtdx);
			/* ionospheric delay (m) */
			if (opt->ionoopt==IONOOPT_EST) {
				vl[nv]-=getIon(x,&HL[nv*nx],ix,UII(i+1,sat,opt),carr,freq);
			}
			/* carrier-phase ambiguity and bias */
			if (carr) {
				vl[nv]-=getAmb(x,&HL[nv*nx],ix,UIAMB(i+1,sat,idx+1,opt),
					rcv[i].datum[sat-1],freq,opt->ionoopt);
				if (uuar) {
					vl[nv]-=getBias(opt,x,&HL[nv*nx],ix,i+1,sat,idx+1,freq);
				}
			}
				
			/* variance */
			var[nv]=varerr(sys,ssat->azel[1],idx,carr,opt)+vart;

			/* record residuals */
			if (!carr) ssat->resp[idx]=vl[nv];	/* code */
			else ssat->resc[idx]=vl[nv];		/* phase */

			/* record large post-fit residuals */
			if (post&&fabs(vl[nv])>sqrt(var[nv])*THRES_REJECT_RATIO) {
				rcvi[ne]=i; obsi[ne]=j; frqi[ne]=k; ve[ne++]=vl[nv];
			}

			if(carr) ssat->vsat[idx]=1;
			nv++;
		}
	}
	/* reject satellite with large and max post-fit residual */
	if (post&&ne>0) {
		postRejSat(uck->time,post,rcv,rcvi,obsi,frqi,ve,ne);
		stat=0;
	}
	/* only generate variance for filter */
	if (!post||flag==LST_EPOCH) {
		if (!post&&flag==2) flag=MID_EPOCH;
		nv=datumConstraint_u(uck,rcv,n,ix,vl,HL,DL,var,nv,flag);
		*nc=uckConstraint(uck,rcv,n,x,P,ix,vp,HP,DP,flag);
	}

	free(rcvi);free(obsi);free(frqi);free(var);free(ve);
	return post?stat:nv;
}