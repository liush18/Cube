
#include "dckp.h"

#define MAX_ITER        8              /* max number of iterations */

/* number of estimated states ------------------------------------------------*/
extern int userPPPnx_dp(const prcopt_t *opt)
{
    return DPNX(opt);
}
/* not estimate datum parameters -----------------------------------------------*/
static void rejEstAmbDatum(const rtk_t *rtk, const int *datum, int *ix)
{
    const prcopt_t *opt=&rtk->opt;
    int i,k,sat;

    if (datum[MAXSAT]) return;

    /* get datum sat */
    for (i=sat=0;i<MAXSAT;i++) if (datum[i]) {
        sat=i+1; break;
    }
    if (!sat) return;

    for (i=0;i<DPNF(opt);i++) {
        k=DPIB(sat,i+1,opt);
        ix[k]=0;
    }
}
/* precise point positioning -------------------------------------------------*/
extern void userPPPos_dp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    const prcopt_t *opt=&rtk->opt;
    double *xp,*Pp,*vl,*HL,*DL,*vp,*HP,*DP;
    double *rs,*dts,*vars,*varc,*azel,dr[3]={0};
    const int nx=rtk->nx;
    int svh[MAXOBS],exc[MAXOBS]={0},datum[MAXSAT+1],stat=SOLQ_SINGLE;
    int i,j,nv,nc,info,*ix;

    trace(3,"userPPPos_dp: time=%s nx=%d n=%d\n",rtk->cTime,nx,n);

    rs=mat(6,n);dts=mat(3,n);vars=mat(1,n);varc=zeros(3,n);azel=zeros(2,n);

    for (i=0;i<MAXSAT;i++) for (j=0;j<opt->nf;j++) rtk->ssat[i].fix[j]=0;

    userUdStates_dp(rtk,obs,n,nav);

    /* satellite positions and clocks */
    sat_poss(obs[0].time,"ppp",obs,n,nav,rs,dts,vars,varc,svh);

    /* exclude measurements of eclipsing satellite (block IIA) */
    testeclipse(obs,n,nav,rs);

    /* earth tides correction */
    tidedisp(gpst2utc(obs[0].time),rtk->x,7,&nav->erp,opt->odisp[0],dr);

    nv=DPCK(opt)*n;nc=DPNX(opt);
    ix=imat(nx,1);xp=mat(nx,1);Pp=mat(nx,nx);
    vl=mat(nv,1);HL=mat(nx,nv);DL=mat(nv,nv);
    vp=mat(nc,1);HP=mat(nx,nc);DP=mat(nc,nc);

    for (i=0;i<MAX_ITER;i++) {

        matcpy(xp,rtk->x,nx,1);
        matcpy(Pp,rtk->P,nx,nx);
        imatcpy(datum,rtk->datum,MAXSAT+1,1);

        /* prefit residuals */
        if (!(nv=userPPPRes_dp(0,obs,n,rs,dts,vars,varc,svh,dr,exc,nav,xp,Pp,
            rtk,vl,HL,DL,vp,HP,DP,azel,ix,&nc))) {
            trace(2,"%s ppp (%d) no valid obs data\n",rtk->cTime,i+1);
            break;
        }

        /* set ambiguity datum for decoupled model */
        userSetAmbDatum(rtk,obs,n,exc,azel,datum);
        /* not estimate datum parameters */
        rejEstAmbDatum(rtk,datum,ix);

        /* measurement update of ekf states */
        if ((info=lsqBlock(xp,Pp,ix,vl,HL,DL,vp,HP,DP,nx,nv,nc))) {
            trace(2,"%s ppp (%d) filter error info=%d\n",rtk->cTime,i+1,info);
            break;
        }

        /* get increment */
        for (j=0;j<nx;j++) if (ix[j]) xp[j]+=rtk->x[j];

        /* postfit residuals */
        if (userPPPRes_dp(i+1,obs,n,rs,dts,vars,varc,svh,dr,exc,nav,xp,Pp,rtk,
            vl,HL,DL,vp,HP,DP,azel,ix,&nc)) {
            matcpy(rtk->x,xp,nx,1);
            matcpy(rtk->P,Pp,nx,nx);
            imatcpy(rtk->datum,datum,MAXSAT+1,1);
            stat=SOLQ_PPP;
            break;
        }
    }
    if (i>=MAX_ITER) {
        trace(2,"%s ppp (%d) iteration overflows\n",rtk->cTime,i);
    }
    if (stat==SOLQ_PPP) {
        if (userAmbResol_dp(rtk,obs,n,exc,azel,xp,Pp)&&
            userPPPRes_dp(MAX_ITER+1,obs,n,rs,dts,vars,varc,svh,dr,exc,nav,xp,Pp,
                rtk,vl,HL,DL,vp,HP,DP,azel,ix,&nc)) {
            matcpy(rtk->xa,xp,nx,1);
            matcpy(rtk->Pa,Pp,nx,nx);
#if 0
            /* hold ambiguity */
            matcpy(rtk->x,xp,nx,1);
            matcpy(rtk->P,Pp,nx,nx);
#endif
            stat=SOLQ_FIX;
        }

        /* update solution status */
        userUdSolution(rtk,obs,n,stat);
        logRclk_dp(rtk);
        logZtd_dp(rtk);
        logAmb_dp(rtk,obs,n);
        logRes_dp(rtk,obs,n);
    }
    free(rs);free(dts);free(vars);free(varc);free(azel);
    free(xp);free(Pp);free(ix);
    free(vl);free(HL);free(DL);
    free(vp);free(HP);free(DP);
}
