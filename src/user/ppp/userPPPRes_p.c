

#include "ppp.h"

#define SQR(x)              ((x)*(x))
#define ROUND(x)			(int)floor((x)+0.5)
#define THRES_REJECT        4.0             /* reject threshold of posfit-res (sigma) */
#define EFACT_GPS_L5		10.0            /* error factor of GPS/QZS L5 */

//#define A1		 SQR(FREQ1)/(SQR(FREQ1)-SQR(FREQ2))
//#define B1		-SQR(FREQ2)/(SQR(FREQ1)-SQR(FREQ2))
//#define A2		 FREQ1/(FREQ1-FREQ2)
//#define B2		-FREQ2/(FREQ1-FREQ2)
//#define A3		 FREQ1/(FREQ1+FREQ2)
//#define B3		 FREQ2/(FREQ1+FREQ2)

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
/* ppp increment constraint ----------------------------------------------------*/
static int pppConstraint(const rtk_t *rtk, const obsd_t *obs, int n, 
    const int *ix, const double *x, const double *P, double *vp, double *HP, 
    double *DP)
{
    const prcopt_t *opt=&rtk->opt;
    const int nx=rtk->nx;
    int i,j,k,l,sat,nc=0,*ind=imat(PNX(opt),1);

    /* position constraint */
    for (k=0;k<3;k++) {
        if (ix[k]&&x[k]!=0.0&&P[k*(nx+1)]!=0.0) {
            for (l=0;l<nx;l++) HP[l+nc*nx]=k==l?1.0:0.0;
            vp[nc]=0.0; ind[nc++]=k;
        }
    }

    /* troposphere constraint */
    if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
        j=PIT(opt);
        for (k=j;k<j+PNT(opt);k++) {
            if (ix[k]&&x[k]!=0.0&&P[k*(nx+1)]!=0.0) {
                for (l=0;l<nx;l++) HP[l+nc*nx]=k==l?1.0:0.0;
                vp[nc]=0.0; ind[nc++]=k;
            }
        }
    }

    /* ionosphere constraint, estimated as white noise if prn<0.0 */
    if (opt->ionoopt==IONOOPT_EST&&opt->prn[2]>=0.0) {
        for (i=0;i<MAXSAT;i++) {
            k=PII(i+1,opt);
            if (ix[k]&&x[k]!=0.0&&P[k*(nx+1)]!=0.0) {
                for (l=0;l<nx;l++) HP[l+nc*nx]=k==l?1.0:0.0;
                vp[nc]=0.0; ind[nc++]=k;
            }
        }
    }

    /* receiver phase bias constraint */
    if (UCK(opt)) {
        for (i=0;i<PNF(&rtk->opt);i++) {
            k=PID(i+1,opt);
            if (ix[k]&&x[k]!=0.0&&P[k*(nx+1)]!=0.0) {
                for (l=0;l<nx;l++) HP[l+nc*nx]=k==l?1.0:0.0;
                vp[nc]=0.0; ind[nc++]=k;
            }
        }
    }

    /* ambiguity constraint */
    for (i=0;i<n;i++) {
        sat=obs[i].sat;
        for (j=0;j<PNF(opt);j++) {
            k=PIB(sat,j+1,opt);
            if (ix[k]&&x[k]!=0.0&&P[k*(nx+1)]!=0.0) {
                for (l=0;l<nx;l++) HP[l+nc*nx]=k==l?1.0:0.0;
                vp[nc]=0.0; ind[nc++]=k;
            } 
        }
    }

    /* variance */
    for (i=0;i<nc;i++) {
        k=ind[i];
        for (j=i;j<nc;j++) {
            l=ind[j];
            DP[i+j*nc]=DP[j+i*nc]=P[k+l*nx];
        }
    }

    free(ind);
    return nc;
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
/* get phase bias ------------------------------------------------------------*/
static double getBias_p(int idx, const prcopt_t *opt, const double *xp, 
    double *HL, int *ix, double freq, double *spb)
{
    const double lam=CLIGHT/freq;
    int k=PID(idx+1,opt);
    HL[k]=lam; ix[k]|=1;
    /* rcv phase bias - sat phase bias */
    return xp[k]*lam-spb[idx]*lam;
}
/* get ambiguity -------------------------------------------------------------*/
static double getAmb_p(const rtk_t *rtk, int idx, const prcopt_t *opt, 
    const double *xp, double *HL, int *ix, int sat, double freq, double *varb)
{
    double bias,lam_NL=1.0,lam=CLIGHT/freq,f1=FREQ1,f2=FREQ2,wl=0.0;
    int k=PIB(sat,idx+1,opt);

    //if (opt->modear>ARMODE_OFF) {
    //    lam_NL=CLIGHT/(FREQ1+FREQ2);
    //    wl=rtk->xw[sat-1];
    //    *varb=rtk->Pw[(sat-1)*(rtk->nw+1)];
    //}

    if (opt->ionoopt==IONOOPT_IFLC) {
        /* conventional ppp model, unit (m or cycle of cnes) */
        bias=lam_NL*xp[k]+CLIGHT*(f2/(SQR(f1)-SQR(f2)))*wl;
        HL[k]=lam_NL; ix[k]|=1;
    }
    else { /* estimation ionosphere, unit (cycle) */
        bias=xp[k]*lam;
        HL[k]=lam; ix[k]|=1;
    }
    return bias;
}
/* get ionosphere delay ------------------------------------------------------*/
static double getIon_p(const double *x, double *HL, int *ix, int k, int carr, 
    double freq)
{
    double C=SQR(FREQ1/freq)*(carr?-1.0:1.0);
    HL[k]=C; ix[k]|=1;
    return C*x[k];
}
/* get troposphere delay -----------------------------------------------------*/
static double getTrop_p(const prcopt_t *opt, double *HL, int *ix, int k, 
    const double *dtdx, double dtrp)
{
    int i;

    for (i=0;i<PNT(opt);i++) { HL[k+i]=dtdx[i]; ix[k+i]|=1; }
    return dtrp;
}
/* get clock state values ----------------------------------------------------*/
static double getClock_p(const prcopt_t *opt, const double *xp, double *HL, 
    int *ix, int sys, const double *dts)
{
    double cdtr,cdts;
    int k;

    /* receiver clock index */
    k=sys==SYS_GLO?2:(sys==SYS_GAL?3:(sys==SYS_CMP?4:1));
    k=PIC(k,opt);
    /* satellite and receiver clock (ns) */
    cdts=CLIGHT*dts[0];
    cdtr=xp[k]*CLIGHT*1E-9;
    HL[k]=CLIGHT*1E-9; ix[k]|=1;

    return cdtr-cdts;
}
/* phase and code residuals --------------------------------------------------*/
extern int userPPPRes_p(int post, const obsd_t *obs, int n, const double *rs,
    const double *dts, const double *var_rs, const double *var_dt, 
    const int *svh, const double *dr, int *exc, const nav_t *nav, 
    const double *xp, double *Pp, rtk_t *rtk, double *vl, double *HL, 
    double *DL, double *vp, double *HP, double *DP, double *azel, int *ix, int *nc)
{
    prcopt_t *opt=&rtk->opt;
    double r,rr[3],pos[3],e[3],L[NFREQ],P[NFREQ],Lc,Pc,freq;
    double trpx[3]={0},dtdx[3],dtrp;
    double vart,varb,var[MAXOBS*NFREQ*2],ve[MAXOBS*NFREQ*2],vmax=0.0;
    double scb[NFREQ]={0},var_scb[NFREQ]={0},spb[NFREQ]={0},var_spb[NFREQ]={0};
    const int nx=rtk->nx;
    int i,j,k,sat,sys,nv=0,stat=1,idx,carr;
    int ne=0,obsi[MAXOBS*NFREQ*2],frqi[MAXOBS*NFREQ*2],maxobs,maxfrq,rej;

    for (i=0;i<MAXSAT;i++) for (j=0;j<opt->nf;j++) rtk->ssat[i].vsat[j]=0;

    for (i=0;i<3;i++) rr[i]=xp[i]+dr[i];
    ecef2pos(rr,pos);

    /* estimate position */
    for (i=0;i<nx;i++) ix[i]=(i<3)?1:0;

    for (i=0;i<n&&i<MAXOBS;i++) {
        sat=obs[i].sat;

        if ((r=geodist(rs+i*6,rr,e))<=0.0||satazel(pos,e,azel+i*2)<opt->elmin) {
            exc[i]=1;
            continue;
        }
        /* relativity correction */
        relativity(rs+i*6,rr,&r);

        if (!(sys=satsys(sat,NULL))||!rtk->ssat[sat-1].vs||
            satexclude(obs[i].sat,var_rs[i],svh[i],opt)||exc[i]) {
            exc[i]=1;
            continue;
        }

        /* tropospheric model */
        matcpy(trpx,xp+PIT(opt),PNT(opt),1);
        tropModel(obs[i].time,pos,azel+i*2,trpx,dtdx,&dtrp,&vart);

        /* calculate corrected measurements */
        CalCorrAntObs(rtk,obs+i,nav,rs+i*6,rr,azel+i*2,L,P,&Lc,&Pc);

        /* mw widelane combination */
        if (UCK(opt)) {
            if (!interpol_phasebias(opt,obs[0].time,sat,nav,spb,var_spb)) 
                continue;
        }
        if (opt->scb) {
            if (interpol_codebias(opt,obs[0].time,sat,nav,scb,var_scb)<0) 
                continue;
        }

        /* stack phase and code residuals {P3,L3,A4,...} or {P1,L1,P2,L2,...} */
        for (j=0;j<2*PNF(opt);j++) {

            varb=0.0;

            idx=j/2;  /* frequency, from 0 */
            carr=j%2; /* obs type, 0 for pseudorange, 1 for carrier-phase */

            if ((freq=sat2freq(sat,obs[i].code[j/2],nav))==0.0) continue;

            /* observables */
            if (opt->ionoopt==IONOOPT_IFLC) {
                if ((vl[nv]=j==0?Pc:Lc)==0.0) continue;
            }
            else {
                if ((vl[nv]=carr==0?P[idx]:L[idx])==0.0) continue;
            }

            vl[nv]-=r;
            if (opt->scb&&!carr) vl[nv]-=-scb[idx]*CLIGHT*1E-9;

            /* coordinates */
            for (k=0;k<nx;k++) HL[k+nv*nx]=(k<3)?-e[k]:0.0;
            /* receiver clock index */
            vl[nv]-=getClock_p(opt,xp,&HL[nv*nx],ix,sys,&dts[i*2]);
            /* troposphere */
            vl[nv]-=getTrop_p(opt,&HL[nv*nx],ix,PIT(opt),dtdx,dtrp);
            /* ionosphere */
            if (opt->ionoopt==IONOOPT_EST) {
                vl[nv]-=getIon_p(xp,&HL[nv*nx],ix,PII(sat,opt),carr,freq);
            }
            /* ambiguity and phase bias */
            if (carr) {
                vl[nv]-=getAmb_p(rtk,idx,opt,xp,&HL[nv*nx],ix,sat,freq,&varb);
                if (UCK(opt)) {
                    vl[nv]-=getBias_p(idx,opt,xp,&HL[nv*nx],ix,freq,spb);
                }
            }

            var[nv]=varerr(sys,azel[1+i*2],idx,carr,opt)+vart+varb;

            /* record residual */
            if (!carr) rtk->ssat[sat-1].resp[idx]=vl[nv];
            else       rtk->ssat[sat-1].resc[idx]=vl[nv];

            trace(4,"%s (%d) sat=%2d %s%d res=%9.4f sig=%9.4f el=%4.1f\n",
                rtk->cTime,post,sat,carr?"L":"P",idx+1,vl[nv],sqrt(var[nv]),
                azel[1+i*2]*R2D);

            /* record large post-fit residuals */
            if (post&&fabs(vl[nv])>sqrt(var[nv])*THRES_REJECT) {
                obsi[ne]=i; frqi[ne]=j; ve[ne]=vl[nv]; ne++;
            }

            if (carr) rtk->ssat[sat-1].vsat[idx]=1;
            nv++;
        }
    }
    /* reject satellite with large and max post-fit residual */
    if (post&&ne>0) {
        vmax=ve[0]; maxobs=obsi[0]; maxfrq=frqi[0]; rej=0;
        for (j=1;j<ne;j++) {
            if (fabs(vmax)>=fabs(ve[j])) continue;
            vmax=ve[j]; maxobs=obsi[j]; maxfrq=frqi[j]; rej=j;
        }
        sat=obs[maxobs].sat;
        trace(2,"%s outlier (%d) rejected sat=%2d %s%d res=%9.4f el=%4.1f\n",
            rtk->cTime,post,sat,maxfrq%2?"L":"P",maxfrq/2+1,vmax,
            azel[1+maxobs*2]*R2D);
        exc[maxobs]=1; rtk->ssat[sat-1].rejc[0]++; stat=0;
        ve[rej]=0;
    }
    if (!post) {
        /* constraint */
        *nc=pppConstraint(rtk,obs,n,ix,xp,Pp,vp,HP,DP);
        for (i=0;i<nv;i++) for (j=0;j<nv;j++) DL[i+j*nv]=(i==j)?var[i]:0.0;
    }

    return post?stat:nv;
}