
#include "dckp.h"

#define SQR(x)      ((x)*(x))
#define THRES_REJECT 4.0            /* reject threshold of posfit-res (sigma) */

#define A1		 SQR(FREQ1)/(SQR(FREQ1)-SQR(FREQ2))
#define B1		-SQR(FREQ2)/(SQR(FREQ1)-SQR(FREQ2))
#define A2		 FREQ1/(FREQ1-FREQ2)
#define B2		-FREQ2/(FREQ1-FREQ2)
#define A3		 FREQ1/(FREQ1+FREQ2)
#define B3		 FREQ2/(FREQ1+FREQ2)

/* ppp increment constraint ----------------------------------------------------*/
static int pppConstraint(const rtk_t *rtk, const obsd_t *obs, int n, 
    const int *ix, const double *x, const double *P, double *vp, double *HP, 
    double *DP)
{
    const prcopt_t *opt=&rtk->opt;
    const int nx=rtk->nx;
    int i,j,k,l,sat,nc=0,*ind=imat(DPNX(opt),1);

    /* position constraint */
    for (k=0;k<3;k++) {
        if (ix[k]&&x[k]!=0.0&&P[k*(nx+1)]!=0.0) {
            for (l=0;l<nx;l++) HP[l+nc*nx]=k==l?1.0:0.0;
            vp[nc]=0.0; ind[nc++]=k;
        }
    }

    /* troposphere constraint */
    if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
        j=DPIT(opt);
        for (k=j;k<j+DPNT(opt);k++) {
            if (ix[k]&&x[k]!=0.0&&P[k*(nx+1)]!=0.0) {
                for (l=0;l<nx;l++) HP[l+nc*nx]=k==l?1.0:0.0;
                vp[nc]=0.0; ind[nc++]=k;
            }
        }
    }

    /* ambiguity constraint */
    for (i=0;i<n;i++) {
        sat=obs[i].sat;
        for (j=0;j<DPNF(opt);j++) {
            k=DPIB(sat,j+1,opt);
            if (ix[k]&&x[k]!=0.0&&P[k*(nx+1)]!=0.0) {
                for (l=0;l<nx;l++) HP[l+nc*nx]=k==l?1.0:0.0;
                vp[nc]=0.0;
                ind[nc++]=k;
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
/* get ambiguity -------------------------------------------------------------*/
static double getAmb_dp(int j, const prcopt_t *opt, const double *xp, double *HL,
    int *ix, int sat, double freq)
{
    double bias;
    int k;

    if (j==0) return 0.0;
    else if (j==1) { /* carrier-phase equation, unit (cycle) */
        k=DPIB(sat,1,opt); 
        bias=xp[k]*17.0*LAM3;
        HL[k]=17.0*LAM3;
        ix[k]|=1;

        k=DPIB(sat,2,opt);
        bias+=xp[k]*60.0*LAM3;
        HL[k]=60.0*LAM3;
        ix[k]|=1;
    }
    else { /* mw equation, unit (cycle) */
        k=DPIB(sat,2,opt);
        bias=xp[k]*LAM4;
        HL[k]=LAM4;
        ix[k]|=1;
    }

    return bias;
}
/* get troposphere delay -----------------------------------------------------*/
static double getTrop_dp(int j, const prcopt_t *opt, double *HL, int *ix, int k, 
    const double *dtdx, double dtrp)
{
    int i;

    if (j==2) return 0.0;
    else {
        for (i=0;i<DPNT(opt);i++) { HL[k+i]=dtdx[i]; ix[k+i]|=1; }
    }
    return dtrp;
}
/* get clock state values ----------------------------------------------------*/
static double getClock_dp(int j, const prcopt_t *opt, const double *xp, 
    double *HL, int *ix, int sys, const double *dts)
{
    double cdtr,cdts;
    int k;

    /* receiver clock index */
    k=sys==SYS_GLO?2:(sys==SYS_GAL?3:(sys==SYS_CMP?4:1));
    k=DPIC(k,j+1,opt);
    ix[k]|=1;
    /* satellite and receiver clock (ns) */
    cdts=j<2?CLIGHT*dts[j]:dts[j]*LAM4;
    cdtr=xp[k]*(j<2?CLIGHT*1E-9:LAM4);
    HL[k]=j<2?CLIGHT*1E-9:LAM4;

    return cdtr-cdts;
}
/* phase and code residuals --------------------------------------------------*/
extern int userPPPRes_dp(int post, const obsd_t *obs, int n, const double *rs,
    const double *dts, const double *var_rs, const double *var_dt, 
    const int *svh, const double *dr, int *exc, const nav_t *nav, 
    const double *xp, double *Pp, rtk_t *rtk, double *vl, double *HL, 
    double *DL, double *vp, double *HP, double *DP, double *azel, int *ix, 
    int *nc)
{
    prcopt_t *opt=&rtk->opt;
    double r,rr[3],pos[3],e[3],L[NFREQ],P[NFREQ],Lc,Pc,A4;
    double trpx[3]={0},dtdx[3],dtrp,vart,tvar,varl,varp,*ve,vmax=0.0,*R;
    char obstype[4];
    const int nx=rtk->nx;
    int i,j,k,sat,sys,stat=1,ne=0,*obsi,*equi,maxobs,maxequ,rej;
    int nv=0,tnv=DPCK(opt)*n;

    for (i=0;i<MAXSAT;i++) for (j=0;j<opt->nf;j++) rtk->ssat[i].vsat[j]=0;

    for (i=0;i<3;i++) rr[i]=xp[i]+dr[i];
    ecef2pos(rr,pos);

    R=zeros(tnv,tnv);
    obsi=izeros(tnv,1);equi=izeros(tnv,1);ve=zeros(tnv,1);

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
        matcpy(trpx,xp+DPIT(opt),DPNT(opt),1);
        tropModel(obs[i].time,pos,azel+i*2,trpx,dtdx,&dtrp,&vart);

        /* calculate corrected measurements */
        CalCorrAntObs(rtk,obs+i,nav,rs+i*6,rr,azel+i*2,L,P,&Lc,&Pc);

        /* mw widelane combination */
        combWideLane(L,P,&A4);

        /* variance of phase and pseudorange */
        varl=SQR(opt->err[1])+SQR(opt->err[2]/sin(azel[i*2+1]));
        varp=SQR(opt->eratio[0])*varl;

        /* stack phase and code residuals {P3,L3,A4,...} */
        for (j=0;j<DPCK(opt);j++) {

            //if ((freq=sat2freq(sat,obs[i].code[j/2],nav))==0.0) continue;

            /* observables */
            if ((vl[nv]=j==0?Pc:(j==1?Lc:A4))==0.0) continue;
            if (j!=2) vl[nv]-=r;

            /* coordinates */
            for (k=0;k<nx;k++) HL[k+nv*nx]=(k<3&&j<2)?-e[k]:0.0;
            /* receiver clock index */
            vl[nv]-=getClock_dp(j,opt,xp,&HL[nv*nx],ix,sys,&dts[i*DPCK(opt)]);
            /* troposphere */
            vl[nv]-=getTrop_dp(j,opt,&HL[nv*nx],ix,DPIT(opt),dtdx,dtrp);
            /* ambiguity */
            vl[nv]-=getAmb_dp(j,opt,xp,&HL[nv*nx],ix,sat,0);

            /* variance, only DCK model has var_dt */
            if (j==0) {
                tvar=SQR(A1)*varp+SQR(B1)*varp+vart+var_rs[i]+var_dt[3*i+j];
            }
            else if (j==1) {
                tvar=SQR(A1)*varl+SQR(B1)*varl+vart+var_rs[i]+var_dt[3*i+j];
            }
            else {
                tvar=SQR(A2)*varl+SQR(B2)*varl+SQR(A3)*varp+SQR(B3)*varp+
                    var_dt[3*i+j];
            }

            R[nv+nv*tnv]=tvar;
            if (j==2) {
                /* cov(P3,A4) */
                R[nv+(nv-2)*tnv]=R[(nv-2)+nv*tnv]=-A1*A3*varp-B1*B3*varp; 
                /* cov(L3,A4) */
                R[nv+(nv-1)*tnv]=R[(nv-1)+nv*tnv]= A1*A2*varl+B1*B2*varl; 
            }

            /* record residual */
            if (j==0)       rtk->ssat[sat-1].resp[0]=vl[nv];
            else if (j==1)  rtk->ssat[sat-1].resc[0]=vl[nv];
            else            rtk->ssat[sat-1].resw[0]=vl[nv];

            sprintf(obstype,"%s",j==0?"P3":(j==1?"L3":"A4"));
            trace(4,"%s (%d) sat=%2d %s res=%9.4f sig=%9.4f el=%4.1f\n",rtk->cTime,post,
                sat,obstype,vl[nv],sqrt(R[nv+nv*tnv]),azel[1+i*2]*R2D);

            /* record large post-fit residuals */
            if (post&&fabs(vl[nv])>sqrt(R[nv+nv*tnv])*THRES_REJECT) {
                obsi[ne]=i; equi[ne]=j; ve[ne]=vl[nv]; ne++;
            }
            if (j==2) rtk->ssat[sat-1].vsat[0]=1;

            nv++;
        }
    }

    /* reject satellite with large and max post-fit residual */
    if (post&&ne>0) {
        vmax=ve[0]; maxobs=obsi[0]; maxequ=equi[0]; rej=0;
        for (j=1;j<ne;j++) {
            if (fabs(vmax)>=fabs(ve[j])) continue;
            vmax=ve[j]; maxobs=obsi[j]; maxequ=equi[j]; rej=j;
        }

        sprintf(obstype,"%s",maxequ==0?"P3":(maxequ==1?"L3":"A4"));

        sat=obs[maxobs].sat;
        trace(2,"%s outlier (%d) rejected sat=%2d %s res=%9.4f el=%4.1f\n",
            rtk->cTime,post,sat,obstype,vmax,azel[1+maxobs*2]*R2D);
        exc[maxobs]=1; rtk->ssat[sat-1].rejc[0]++; stat=0;
        ve[rej]=0;
    }
    if (!post) {
        /* constraint */
        *nc=pppConstraint(rtk,obs,n,ix,xp,Pp,vp,HP,DP);

        for (i=0;i<nv;i++) {
            for (j=0;j<nv;j++) DL[i+j*nv]=R[i+j*tnv];
        }
    }

    free(R);free(obsi);free(equi);free(ve);
    return post?stat:nv;
}