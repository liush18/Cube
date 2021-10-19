

#include "ppp.h"

#define SQR(x)          ((x)*(x))
#define ROUND(x)        (int)floor((x)+0.5)
#define MAX_SAT_OUTC    5               /* obs outage count to reset bias */
#define MAX_ITER        8               /* max number of iterations */
#define THRES_REJECT    4.0             /* reject threshold of posfit-res (sigma) */
#define MIN_ARC_GAP     300.0           /* min arc gap (s) */
/* not estimate datum parameters -----------------------------------------------*/
static void rejMWDatumEst(const rtk_t *rtk, const int *datum, int *ix)
{
    int i,sat;

    if (datum[MAXSAT]) return;

    /* get datum sat */
    for (i=sat=0;i<MAXSAT;i++) if (datum[i]) {
        sat=i+1; break;
    }
    if (!sat) return;

    ix[sat-1]=0;
}
/* wide-lane ambiguity constraint --------------------------------------------*/
static int MWConstraint(const rtk_t *rtk, const obsd_t *obs, int n, 
    const int *ix, const double *x, const double *P, double *vp, double *HP, 
    double *DP)
{
    const int nw=rtk->nw;
    int i,j,k,l,nc,*ind=imat(MAXSAT+1,1);

    for (i=nc=0;i<MAXSAT+1;i++) {
        if (ix[i]&&x[i]!=0.0&&P[i*(nw+1)]!=0.0) {
            for (j=0;j<nw;j++) HP[j+nc*nw]=i==j?1.0:0.0;
            vp[nc]=0.0; ind[nc++]=i;
        }
    }

    /* variance */
    for (i=0;i<nc;i++) {
        k=ind[i];
        for (j=i;j<nc;j++) {
            l=ind[j];
            DP[i+j*nc]=DP[j+i*nc]=P[k+l*nw];
        }
    }

    free(ind);
    return nc;
}
/* measurement error variance ------------------------------------------------*/
static double varerr_mw(double el, const prcopt_t *opt)
{
    double fact=3.0,sinel=sin(el),varp,varl;
    double f1=FREQ1,f2=FREQ2;

    varl=SQR(fact*opt->err[1])+SQR(fact*opt->err[2]/sinel);
    varp=varl*opt->eratio[0];

    return (SQR(f1/(f1-f2))+SQR(f2/(f1-f2)))*varl+
        (SQR(f1/(f1+f2))+SQR(f2/(f1+f2)))*varp;
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
/* mw residuals --------------------------------------------------------------*/
static int MWRes(int post, rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
    const double *rs, const double *var_rs, const int *svh, const double *dr, 
    const double *xw, const double *Pw,  double *vl, double *HL,  double *DL, 
    double *vp, double *HP, double *DP, double *azel, int *exc, int *ix, int *nc)
{
    const prcopt_t *opt=&rtk->opt;
    ambc_t *amb;
    double f1=FREQ1,f2=FREQ2,wl=CLIGHT/(f1-f2),L[NFREQ],P[NFREQ],Lc,Pc,LC1,var1;
    double rr[3],pos[3],e[3],var[MAXOBS],ve[MAXOBS],vmax;
    const int nw=rtk->nw;
    int i,j,nv,ne,sys,sat,obsi[MAXOBS],maxobs,rej,stat=1;

    double dant[NFREQ]={0};

    for (i=0;i<3;i++) rr[i]=rtk->x[i]+dr[i];
    for (i=0;i<nw;i++) ix[i]=0;
    ecef2pos(rr,pos);

    for (i=nv=ne=0;i<n&&i<MAXOBS;i++) {
        sat=obs[i].sat;
        if (geodist(rs+i*6,rr,e)<=0.0||satazel(pos,e,azel+i*2)<opt->elmin) {
            exc[i]=1;
            continue;
        }

        if (!(sys=satsys(sat,NULL))||!rtk->ssat[sat-1].vs||
            satexclude(obs[i].sat,var_rs[i],svh[i],opt)||exc[i]) {
            exc[i]=1;
            continue;
        }

        for (j=0;j<nw;j++) HL[j+nv*nw]=0.0;

        amb=rtk->ambc+sat-1;
        if (!post) {
            /* calculate corrected measurements */
            //CalCorrAntObs(rtk,obs+i,nav,rs+i*6,rr,azel+i*2,L,P,&Lc,&Pc);
            corr_meas(obs+i,nav,azel,opt,dant,dant,0.0,L,P,&Lc,&Pc);
            LC1=L_LC(1,-1, 0,L)-P_LC(1,1,0,P);
            var1=varerr_mw(azel[1+i*2],&rtk->opt);

            /* reset */
            if (rtk->ssat[sat-1].slip[0]||rtk->ssat[sat-1].slip[1]||amb->n[0]==0
                ||fabs(timediff(amb->epoch[0],obs[0].time))>MIN_ARC_GAP) {

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

        vl[nv]=amb->LC[0];
        var[nv]=amb->LCv[0];

        //if (post) trace(2,"%s sat=G%02d wl=%8.3f\n",rtk->cTime,sat,vl[nv]);

        vl[nv]-=wl*(xw[sat-1]+xw[MAXSAT]-nav->wlbias[sat-1]);
        HL[sat-1+nv*nw]=wl;
        ix[sat-1]|=1;
        HL[MAXSAT+nv*nw]=wl;
        ix[MAXSAT]|=1;

        rtk->ssat[sat-1].resw[0]=vl[nv];
        if (post&&fabs(vl[nv])>sqrt(var[nv])*THRES_REJECT) {
            obsi[ne]=i; ve[ne]=vl[nv]; ne++;
        }

        nv++;
    }
    if (post&&ne>0) {
        vmax=ve[0]; maxobs=obsi[0]; rej=0;
        for (i=1;i<ne;i++) {
            if (fabs(vmax)>=fabs(ve[i])) continue;
            vmax=ve[i]; maxobs=obsi[i], rej=i;
        }
        sat=obs[maxobs].sat;
        trace(2, "%s outlier (%d) wl rejected sat=%02d res=%9.4f el=%4.1f\n",
            rtk->cTime,post,sat,vmax,azel[1+maxobs*2]*R2D);
        exc[maxobs]=1; stat=0; ve[rej]=0;
    }
    if (!post) {
        *nc=MWConstraint(rtk,obs,n,ix,xw,Pw,vp,HP,DP);
        for (i=0;i<nv;i++) for (j=0;j<nv;j++) DL[i+j*nv]=(i==j)?var[i]:0.0;
    }

    return post?stat:nv;
}
static void initMWx(rtk_t *rtk, double xi, double var, int i)
{
    int j;
    rtk->xw[i]=xi;
    for (j=0;j<rtk->nw;j++) {
        rtk->Pw[i+j*rtk->nw]=rtk->Pw[j+i*rtk->nw]=i==j?var:0.0;
    }
}
/* temporal update of mw ambiguity and bias ----------------------------------*/
static void udMW(rtk_t *rtk, const obsd_t *obs, int n)
{
    prcopt_t *opt=&rtk->opt;
    int i,sat,slip,clk_jump=0;

    clk_jump=ROUND(time2gpst(obs[0].time,NULL)*10)%864000==0;

    /* reset receiver mw bias */
    if (clk_jump) initMWx(rtk,0.0,0.0,MAXSAT);
    for (i=0;i<MAXSAT;i++) {
        if (++rtk->ssat[i].outc[0]>MAX_SAT_OUTC||
            opt->modear==ARMODE_INST||clk_jump) {
            initMWx(rtk,0.0,0.0,i);
        }
    }

    for (i=0;i<n&&i<MAXOBS;i++) {

        sat=obs[i].sat;
        slip=rtk->ssat[sat-1].slip[0]||rtk->ssat[sat-1].slip[1];

        /* x!=0.0 and not slip */
        if (rtk->xw[sat-1]!=0.0&&!slip) continue; 
        initMWx(rtk,0.0,0.0,sat-1);
    }
}
/* estimate and fix mw ambiguity ---------------------------------------------*/
extern void userMWEst_cnes(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav, 
    const double *rs, const double *var_rs, const int *svh, const double *dr, 
    int *exc)
{
    double *xw,*Pw,*vl,*HL,*DL,*vp,*HP,*DP,*azel;
    const int nw=rtk->nw;
	int i,j,nv,nc,*ix,info,stat=0,datum[MAXSAT+1];

	trace(3,"userMWEst_cnes: time=%s n=%d\n",rtk->cTime,n);

    azel=zeros(2,n);

    udMW(rtk,obs,n);

    nv=n;nc=n+1;
    ix=imat(nw,1);xw=mat(nw,1);Pw=mat(nw,nw);
    vl=mat(nv,1);HL=mat(nw,nv);DL=mat(nv,nv);
    vp=mat(nc,1);HP=mat(nw,nc);DP=mat(nc,nc);
    for (i=0;i<MAX_ITER;i++) {

        matcpy(xw,rtk->xw,nw,1);
        matcpy(Pw,rtk->Pw,nw,nw);
        imatcpy(datum,rtk->datum,MAXSAT+1,1);

        if (!(nv=MWRes(0,rtk,obs,n,nav,rs,var_rs,svh,dr,xw,Pw,vl,HL,DL,vp,HP,DP,
            azel,exc,ix,&nc))) {
            trace(2,"%s userMWEst_cnes (%d) no valid data\n",rtk->cTime,i+1);
        }

        userSetAmbDatum(rtk,obs,n,exc,azel,datum);
        rejMWDatumEst(rtk,datum,ix);

        if ((info=lsqBlock(xw,Pw,ix,vl,HL,DL,vp,HP,DP,nw,nv,nc))) {
            trace(2,"%s ppp (%d) filter error info=%d\n",rtk->cTime,i+1,info);
            break;
        }

        /* get increment */
        for (j=0;j<nw;j++) if (ix[j]) xw[j]+=rtk->xw[j];

        if (nv=MWRes(i+1,rtk,obs,n,nav,rs,var_rs,svh,dr,xw,Pw,vl,HL,DL,vp,HP,DP,
            azel,exc,ix,&nc)) {
            matcpy(rtk->xw,xw,nw,1);
            matcpy(rtk->Pw,Pw,nw,nw);
            imatcpy(rtk->datum,datum,MAXSAT+1,1);
            stat=1;
            break;
        }
    }
    if (stat) {
        if (userMWAmbResol_cnes(rtk,obs,n,exc,xw,Pw)&&
            MWRes(MAX_ITER,rtk,obs,n,nav,rs,var_rs,svh,dr,xw,Pw,vl,HL,DL,vp,HP,
                DP,azel,exc,ix,&nc)) {
            matcpy(rtk->xw,xw,nw,1);
            matcpy(rtk->Pw,Pw,nw,nw);
        }
        //trace(2,"%s\n",rtk->cTime);
        //for (i=0;i<nw;i++) if (ix[i]) {
        //    trace(2,"sat=G%02d wl=%10.3f\n",i+1,xw[i]);
        //}
    }

    free(azel);
    free(xw);free(Pw);free(ix);
    free(vl);free(HL);free(DL);
    free(vp);free(HP);free(DP);
}