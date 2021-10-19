
#include "dck.h"

#define SQR(x)	 (x*x)
#define A1		 SQR(FREQ1)/(SQR(FREQ1)-SQR(FREQ2))
#define B1		-SQR(FREQ2)/(SQR(FREQ1)-SQR(FREQ2))
#define A2		 FREQ1/(FREQ1-FREQ2)
#define B2		-FREQ2/(FREQ1-FREQ2)
#define A3		 FREQ1/(FREQ1+FREQ2)
#define B3		 FREQ2/(FREQ1+FREQ2)
#define MIN_VALUE_ZERO		1E-10
#define CONST_SAT_ZMC		1E-3	/* constraint to ZMC */
#define THRES_REJECT_RATIO	4.0		/* threshold ratio of variance to reject */

/* constraint by predicted values --------------------------------------------*/
static int dckConstraint(const net_t *dck, const rcv_t *rcv, int n, 
	const double *x, const double *P, const int *ix, double *vp, double *HP, 
	double *DP, int flag)
{
	const prcopt_t *opt=&dck->opt;
	const int nx=dck->nx;
	int i,j,k,l,f,sat,nc=0,*id;

	id=imat(nx,1);
	openRecord_d(flag);
	record_d(TYPE_TIME,dck->time,"",0,0,0,0,0.0,flag);

	for (i=0;i<n&&i<MAXRCV;i++) {

		if (rcv[i].nobs==0) continue;

		/* record tropospheric delay constraint info */
		k=DIT(i+1,opt);
		for (j=0;j<DNT(opt);j++) {
			if (ix[k+j]&&x[k+j]!=0.0&&P[(k+j)*(nx+1)]!=0.0) {
				record_d(TYPE_TRP,"",rcv[i].sta.name,0,0,j,0,x[k+j],flag);
				for (l=0;l<nx;l++) HP[l+nc*nx]=(k+j)==l?1.0:0.0;
				vp[nc]=0.0;
				id[nc++]=k+j;
			}
		}

		/* record ambiguity constraint info */
		for (j=0;j<rcv[i].nobs&&j<MAXOBS;j++) {

			sat=rcv[i].obs[j].sat;
			if (rcv[i].datum[sat-1]<OBS_USE) continue;

			/* record ambiguity datum */
			if (rcv[i].datum[sat-1]>OBS_USE) {
				k=DIAMB(i+1,sat,1,opt); l=DIAMB(i+1,sat,2,opt);
				record_d(TYPE_AMB,"",rcv[i].sta.name,sat,1,0,
					rcv[i].datum[sat-1],x[k],flag);
				record_d(TYPE_AMB,"",rcv[i].sta.name,sat,2,0,
					rcv[i].datum[sat-1],x[l],flag);
				continue;
			}

			k=DIAMB(i+1,sat,1,opt); l=DIAMB(i+1,sat,2,opt);
			if (ix[k]&&x[k]!=0.0&&P[k+k*nx]!=0.0&&
				ix[l]&&x[l]!=0.0&&P[l+l*nx]!=0.0) {

				record_d(TYPE_AMB,"",rcv[i].sta.name,sat,1,0,
					rcv[i].datum[sat-1],x[k],flag);
				for (f=0;f<nx;f++) HP[f+nc*nx]=k==f?1.0:0.0;
				id[nc]=k;
				vp[nc++]=0.0;

				record_d(TYPE_AMB,"",rcv[i].sta.name,sat,2,0,
					rcv[i].datum[sat-1],x[l],flag);
				for (f=0;f<nx;f++) HP[f+nc*nx]=l==f?1.0:0.0;
				id[nc]=l;
				vp[nc++]=0.0;
			}
		}
	}
	/* trop and ambiguities constraint */
	for (i=0;i<nc;i++) {
		k=id[i];
		for (j=i;j<nc;j++) {
			l=id[j];
			DP[i+j*nc]=DP[j+i*nc]=P[k+l*nx];
		}
	}

	recordVar_d(DP,nc,flag);
	closeRecord_d(flag);
	free(id);
	return nc;
}
/* datum costraint eliminate rank-deficiency ---------------------------------*/
static int datumConstraint_d(const net_t *dck, const rcv_t *rcv, int n, 
	const int *ix, double *vl, double *HL, double *DL, int nv, int tnv, 
	int epFlag)
{
	const prcopt_t *opt=&dck->opt;
	double *DL_;
	int i,j,k,nx=dck->nx;

	/* zero-mean condition imposed on all satellite clocks */
	if (opt->datumType==DTM_ZMCSAT&&opt->estclk) {
		for (i=0;i<CLK_NUM(opt);i++) {
			for (j=0;j<nx;j++) HL[j+nv*nx]=0.0;
			for (j=0;j<MAXSAT;j++) {
				k=DISC(j+1,i+1,opt);
				if (ix[k]) HL[k+nv*nx]=1.0;
			}
			vl[nv]=0.0;
			DL[nv+nv*tnv]=SQR(CONST_SAT_ZMC);
			nv++;
		}
	}

	/* measurement variance */
	DL_=mat(nv,nv);
	for (i=0;i<nv;i++) for (j=0;j<nv;j++) DL_[i+j*nv]=DL[i+j*tnv];
	matcpy(DL,DL_,nv,nv);
	
	free(DL_);
	return nv;
}
/* reject satellite with large and max post-fit residual ---------------------*/
static void postRejSat_d(int post, const char *stime, rcv_t *rcv, int n,
	const int *rcvi, const int *obsi, const int *equi, double *ve, int ne,
	int *outr, int *outs)
{
	double vmax;
	char id[8];
	int j,sat,maxrcv,maxobs,maxequ,rej;

	vmax=ve[0]; maxrcv=rcvi[0]; maxobs=obsi[0]; maxequ=equi[0]; rej=0;
	for (j=1;j<ne;j++) {
		if (fabs(vmax)>=fabs(ve[j])) continue;
		vmax=ve[j]; maxrcv=rcvi[j]; maxobs=obsi[j]; maxequ=equi[j]; rej=j;
	}
	sat=rcv[maxrcv].obs[maxobs].sat; satno2id(sat,id);
	trace(2,"%s outlier (%2d) rejected rcv(%4s) sat(%s) %s res=%9.4f el=%4.1f\n",
		stime,post,rcv[maxrcv].sta.name,id,maxequ==0?"Pc":(maxequ==1?"Lc":"Ac"),
		vmax,rcv[maxrcv].ssat[sat-1].azel[1]*R2D);
	rcv[maxrcv].datum[sat-1]=OBS_REJ;
	rcv[maxrcv].ssat[sat-1].rejc[0]++;
	rcv[maxrcv].ssat[sat-1].vsat[0]=0;
	ve[rej]=0;
	/* not reject in ar res */
	if (post!=n) {
		outr[maxrcv]++; outs[sat-1]++;
	}
}
/* get ambiguity -------------------------------------------------------------*/
static double getAmb(const prcopt_t *opt, const double *x, double *HL, int *ix, 
	int ircv, int sat, int k, int datum)
{	
	double bias;
	int j;

	/* ambiguity */
	if (k==0) bias=0.0; /* code */
	else if (k==1) {	/* phase */
		j=DIAMB(ircv,sat,1,opt);
		bias=x[j]*LAM3*17.0;
		HL[j]=LAM3*17.0;
		if (datum==OBS_USE) ix[j]|=1;

		j=DIAMB(ircv,sat,2,opt);
		bias+=x[j]*LAM3*60.0;
		HL[j]=LAM3*60.0;
		if (datum==OBS_USE) ix[j]|=1;
	}
	else { /* wide */
		j=DIAMB(ircv,sat,2,opt);
		bias=x[j]*LAM4;
		HL[j]=LAM4;
		if (datum==OBS_USE) ix[j]|=1;
	}

	return bias;
}
/* get troposphere delay -----------------------------------------------------*/
static double getTrop_d(const prcopt_t *opt, double *HL, int *ix, int k, 
	const double *dtdx)
{
	int i;

	for (i=0;i<DNT(opt);i++) {
		HL[k+i]=dtdx[i];
		ix[k+i]|=1;
	}
	return 0.0;
}
/* get clock state values ----------------------------------------------------*/
static double getClock_d(const prcopt_t *opt, const double *x, double *HL, 
	int *ix, int ircv, int sat, int k, double dts)
{
	double cdtr,cdts;
	int j;

	/* receiver clock (ns) */
	if (opt->datumType==DTM_REFRCV&&ircv==opt->ircv) cdtr=0.0;
	else {
		j=DIRC(ircv,k+1,opt);
		cdtr=x[j]*(k<2?CLIGHT*1E-9:LAM4);
		HL[j]=k<2?CLIGHT*1E-9:LAM4;
		ix[j]|=1;
	}

	/* satellite clock (ns) */
	if (!opt->estclk) cdts=-CLIGHT*dts;
	else {
		j=DISC(sat,k+1,opt);
		cdts=-x[j]*(k<2?CLIGHT*1E-9:LAM4);
		HL[j]=k<2?-CLIGHT*1E-9:-LAM4;
		ix[j]|=1;
	}

	return cdtr+cdts;
}
/* decoupled clock model residuals -------------------------------------------*/
extern int dckRes(int post, const net_t *dck, rcv_t *rcv, int n,
	const nav_t *nav, const double *x, const double *P, int *ix, double *vl, 
	double *HL, double *DL, double *vp, double *HP, double *DP, int *nc, 
	int tnv, int flag, int *outr, int *outs)
{
	const prcopt_t *opt=&dck->opt;
	ssat_t *ssat;
	double trpx[3]={0},dtdx[3],dtrp,vart,varp,varl,var_rs,tvar,*ve;
	const int nx=dck->nx;
	int i,j,k,l,sat,sys,stat=1,nv=0,ne=0,*rcvi,*obsi,*equi;

	rcvi=imat(tnv,1); obsi=imat(tnv,1); equi=imat(tnv,1); ve=mat(tnv,1);
	for (i=0;i<nx;i++) ix[i]=0;
	for (i=0;i<tnv*tnv;i++) DL[i]=0.0;
	for (i=0;i<n;i++) for (j=0;j<MAXSAT;j++) rcv[i].ssat[j].vsat[0]=0;

	for (i=0;i<n&&i<MAXRCV;i++) for(j=0;j<rcv[i].nobs&&j<MAXOBS;j++) {

		sat=rcv[i].obs[j].sat;
		ssat=&rcv[i].ssat[sat-1];
		sys=satsys(sat,NULL);

		/* skip rejected satellite */
		if (rcv[i].datum[sat-1]<OBS_USE) continue;

		/* tropospheric delay */
		matcpy(trpx,x+DIT(i+1,opt),DNT(opt),1);
		tropModel(rcv[i].obs[j].time,rcv[i].pos,ssat->azel,trpx,dtdx,&dtrp,&vart);

		var_rs=rcv[i].varrs[j];
		varl=SQR(opt->err[1])+SQR(opt->err[2]/sin(ssat->azel[1]));
		varp=SQR(opt->eratio[0])*varl;
			
		/* pseudorange, carrier phase, mw widelane ... */
		for(k=0;k<CLK_NUM(opt);k++){

			/* obs measurement */
			vl[nv]=k==0?rcv[i].Pc[j]:(k==1?rcv[i].Lc[j]:rcv[i].A4[j]);

			/* initialization design matrix */
			for (l=0;l<nx;l++) HL[l+nv*nx]=0.0;
			/* receiver-satellite distance */
			if (k!=2) vl[nv]-=rcv[i].r[j]+dtrp;
			/* receiver and satellite clock, dts sets to 0.0 */
			vl[nv]-=getClock_d(opt,x,&HL[nv*nx],ix,i+1,sat,k,0.0);
			/* tropospheric delay (m) */
			if (k!=2) vl[nv]-=getTrop_d(opt,&HL[nv*nx],ix,DIT(i+1,opt),dtdx);
			/* ambiguity */
			vl[nv]-=getAmb(opt,x,&HL[nv*nx],ix,i+1,sat,k,rcv[i].datum[sat-1]);

			/* variance */
			if      (k==0) tvar=SQR(A1)*varp+SQR(B1)*varp+var_rs+vart;
			else if (k==1) tvar=SQR(A1)*varl+SQR(B1)*varl+var_rs+vart;
			else           tvar=SQR(A2)*varl+SQR(B2)*varl+
								SQR(A3)*varp+SQR(B3)*varp;
			DL[nv+nv*tnv]=tvar;
			if (k==2) { /* cov(P3,A4) and cov(L3,A4) */
				DL[nv+(nv-2)*tnv]=DL[(nv-2)+nv*tnv]=-A1*A3*varp-B1*B3*varp;
				DL[nv+(nv-1)*tnv]=DL[(nv-1)+nv*tnv]= A1*A2*varl+B1*B2*varl;
			}

			if (k==0) {
				rcv[i].ssat[sat-1].resp[0]=vl[nv];
				rcv[i].ssat[sat-1].resp[1]=tvar;
			}
			else if (k==1) {
				rcv[i].ssat[sat-1].resc[0]=vl[nv];
				rcv[i].ssat[sat-1].resc[1]=tvar;
			}
			else {
				rcv[i].ssat[sat-1].resw[0]=vl[nv];
				rcv[i].ssat[sat-1].resw[1]=tvar;
			}

			/* record large post-fit residuals */
			if (post&&fabs(vl[nv])>sqrt(tvar)*THRES_REJECT_RATIO) {
				rcvi[ne]=i; obsi[ne]=j; equi[ne]=k; ve[ne++]=vl[nv];
			}

			if(k==2) rcv[i].ssat[sat-1].vsat[0]=1;
			nv++;
		}
	}
	/* reject satellite with large and max post-fit residual */
	if (post&&ne>0) {
		postRejSat_d(post,dck->time,rcv,n,rcvi,obsi,equi,ve,ne,outr,outs);
		stat=0;
	}

	if (!post||flag==LST_EPOCH) {
		if (!post&&flag==LST_EPOCH) flag=MID_EPOCH;
		nv=datumConstraint_d(dck,rcv,n,ix,vl,HL,DL,nv,tnv,flag);
		//tracemat(2,DL,nv,nv,10,6);
		*nc=dckConstraint(dck,rcv,n,x,P,ix,vp,HP,DP,flag);
		//tracemat(2,DP,*nc,*nc,10,6);
	}

	free(rcvi); free(obsi); free(equi); free(ve);
	return post?stat:nv;
}