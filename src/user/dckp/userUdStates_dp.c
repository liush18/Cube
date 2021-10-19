

#include "dckp.h"

#define ROUND(x)    (int)floor((x)+0.5)

#define SQR(x)		    (x*x)
#define VAR_POS         SQR(100.0)  /* init variance position (m^2) */
#define VAR_GRA         SQR(0.01)   /* init variance gradient (m^2) */
#define GAP_RESION      120         /* default gap to reset ionos parameters (ep) */
#define MAX_SAT_OUTC    5           /* obs outage count to reset bias */


/* temporal update of position -----------------------------------------------*/
static void udPos(rtk_t *rtk)
{
    int i;

    /* fixed mode */
    if (rtk->opt.mode==PMODE_PPP_FIXED) {
        for (i=0;i<3;i++) userInitX_dp(rtk,rtk->opt.ru[i],1E-8,i);
        return;
    }
    /* initialize position for first epoch */
    if (norm(rtk->x,3)<=0.0) {
        for (i=0;i<3;i++) userInitX_dp(rtk,rtk->sol.rr[i],VAR_POS,i);
    }
    /* static ppp mode */
    if (rtk->opt.mode==PMODE_PPP_STATIC) {
        for (i=0;i<3;i++) {
            rtk->P[i*(1+rtk->nx)]+=rtk->opt.prn[0]*fabs(rtk->tt);
        }
        return;
    }
    /* kinmatic mode */
    for (i=0;i<3;i++) userInitX_dp(rtk,rtk->sol.rr[i],VAR_POS,i);
}
/* temporal update of tropospheric parameters --------------------------------*/
static void udTrop(rtk_t *rtk)
{
    const prcopt_t *opt=&rtk->opt;
    double pos[3],azel[]={0.0,PI/2.0},ztd,var;
    int j,k=DPIT(&rtk->opt),nx=rtk->nx;

    if (rtk->x[k]==0.0) {
        ecef2pos(rtk->sol.rr,pos);
        /* zenith total tropospheric delay (m) */
        ztd=tropMops(rtk->sol.time,pos,azel,&var);
        userInitX_dp(rtk,ztd,var,k);
        if (opt->tropopt==TROPOPT_ESTG) {
            for (j=k+1;j<k+3;j++) userInitX_dp(rtk,1E-6,VAR_GRA,j);
        }
    }
    else {
        rtk->P[k*(nx+1)]+=opt->prn[1]*fabs(rtk->tt);
        if (opt->tropopt==TROPOPT_ESTG) {
            for (j=k+1;j<k+3;j++) {
                rtk->P[j*(nx+1)]+=SQR(0.1)*opt->prn[1]*fabs(rtk->tt);
            }
        }
    }
}
/* temporal update of ambiguity ----------------------------------------------*/
static void udAmb(rtk_t *rtk, const obsd_t *obs, int n)
{
    prcopt_t *opt=&rtk->opt;
    int i,k1,k2,sat,slip,clk_jump=0;

    /* handle day-boundary clock jump, our decoupled clocks are not is 
    discontinuous, so clk_jump is necessary */
    clk_jump=ROUND(time2gpst(obs[0].time,NULL)*10)%864000==0;

    /* reset phase-bias if expire obs outage counter */
    for (i=0;i<MAXSAT;i++) {
        if (++rtk->ssat[i].outc[0]>MAX_SAT_OUTC||opt->modear==ARMODE_INST||
            clk_jump) {
            userInitX_dp(rtk,0.0,0.0,DPIB(i+1,1,opt)); /* narr */
            userInitX_dp(rtk,0.0,0.0,DPIB(i+1,2,opt)); /* wide */
        }
    }

    for (i=0;i<n&&i<MAXOBS;i++) {

        sat=obs[i].sat;
        slip=rtk->ssat[sat-1].slip[0]||rtk->ssat[sat-1].slip[1];

        k1=DPIB(sat,1,opt); k2=DPIB(sat,2,opt);
        if ((rtk->x[k1]!=0.0)&&(rtk->x[k2]!=0.0)&&!slip) {
            rtk->P[k1*(rtk->nx+1)]+=rtk->opt.prn[3]*fabs(rtk->tt);
            rtk->P[k2*(rtk->nx+1)]+=rtk->opt.prn[3]*fabs(rtk->tt);
            continue; 
        }
        userInitX_dp(rtk,0.0,0.0,k1);
        userInitX_dp(rtk,0.0,0.0,k2);
    }
}
/* temporal update of states -------------------------------------------------*/
extern void userUdStates_dp(rtk_t *rtk, const obsd_t *obs, int n, 
    const nav_t *nav)
{
    trace(3,"userUdStates_dp: n=%d\n",n);

    /* temporal update of position */
    udPos(rtk);

    /* temporal update of troposphere */
    if (rtk->opt.tropopt==TROPOPT_EST||rtk->opt.tropopt==TROPOPT_ESTG) {
        udTrop(rtk);
    }

    /* temporal update of ambiguity */
    udAmb(rtk,obs,n);
}