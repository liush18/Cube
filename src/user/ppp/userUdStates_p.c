

#include "ppp.h"

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
    int j,k=PIT(&rtk->opt),nx=rtk->nx;

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
/* temporal update of ionospheric parameters -----------------------------------
* notes : constrains in lsq filter determine the parameter estimation method
*         in udIon, ion can be estimated as white noise or random walk
*         stats-prnion < 0.0 means ion is esitmated as white noise
*         more details in userPPPRes_p
* ----------------------------------------------------------------------------*/
static void udIon(rtk_t *rtk, const obsd_t *obs, int n)
{
    int i,k;

    for (i=0;i<n&&i<MAXOBS;i++) {
        k=PII(obs[i].sat,&rtk->opt);
        if (rtk->x[k]!=0.0&&(int)rtk->ssat[obs[i].sat-1].outc[0]>GAP_RESION) {
            userInitX_dp(rtk,0.0,0.0,k);
        }
        else rtk->P[k*(rtk->nx+1)]+=rtk->opt.prn[2]*fabs(rtk->tt);
    }
}
/* temporal update of ambiguity ----------------------------------------------*/
static void udAmb(rtk_t *rtk, const obsd_t *obs, int n)
{
    prcopt_t *opt=&rtk->opt;
    int i,j,f,sat,slip,clk_jump=0;

    /* handle day-boundary clock jump, our decoupled clocks are not is 
       discontinuous, so clk_jump is necessary */
    clk_jump=ROUND(time2gpst(obs[0].time,NULL)*10)%864000==0;

    for (f=0;f<PNF(opt);f++) {
        /* reset phase-bias if expire obs outage counter */
        for (i=0;i<MAXSAT;i++) {
            if (++rtk->ssat[i].outc[f]>MAX_SAT_OUTC||
                opt->modear==ARMODE_INST||clk_jump) {
                userInitX_dp(rtk,0.0,0.0,PIB(i+1,f+1,opt));
            }
        }

        for (i=0;i<n&&i<MAXOBS;i++) {

            sat=obs[i].sat;
            if (opt->ionoopt==IONOOPT_IFLC) {
                slip=rtk->ssat[sat-1].slip[0]||rtk->ssat[sat-1].slip[1];
            }
            else slip=rtk->ssat[sat-1].slip[f];

            j=PIB(sat,f+1,opt);
            /* x!=0.0 and not slip */
            if (rtk->x[j]!=0.0&&!slip) {
                rtk->P[j*(rtk->nx+1)]+=rtk->opt.prn[3]*fabs(rtk->tt);
                continue; 
            }
            userInitX_dp(rtk,0.0,0.0,j);
        }
    }
}
/* temporal update of receiver phase bias ------------------------------------*/
static void udBias(rtk_t *rtk)
{
    int i,j;

    for (i=0;i<PND(&rtk->opt);i++) {
        j=PID(i+1,&rtk->opt);
        if (rtk->x[j]!=0.0) {
            rtk->P[j*(rtk->nx+1)]+=rtk->opt.prn[4]*fabs(rtk->tt);
        }
    }
}
/* temporal update of states -------------------------------------------------*/
extern void userUdStates_p(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    trace(3,"userUdStates_p: n=%d\n",n);

    /* temporal update of position */
    udPos(rtk);

    /* temporal update of troposphere */
    if (rtk->opt.tropopt==TROPOPT_EST||
        rtk->opt.tropopt==TROPOPT_ESTG) udTrop(rtk);

    /* temporal update of ionosphere */
    if (rtk->opt.ionoopt==IONOOPT_EST) udIon(rtk,obs,n);

    /* temporal update of ambiguity */
    udAmb(rtk,obs,n);

    /* temporal update of phase-bias */
    if (UCK(&rtk->opt)) udBias(rtk);
}