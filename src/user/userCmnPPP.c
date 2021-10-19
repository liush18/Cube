
#include "user.h"

#define MIN_NSAT_SOL    4              /* min satellite number for solution */

/* set ambiguity datum or check datum for ppp --------------------------------*/
extern void userSetAmbDatum(const rtk_t *rtk, const obsd_t *obs, int n, 
    const int *exc, const double *azel, int *datum)
{
    double satel[MAXOBS]={0};
    int i,j,sat,slip;

    if (datum[MAXSAT]) return;

    /* get datum sat */
    for (i=sat=0;i<MAXSAT;i++) if (datum[i]) {
        sat=i+1; break;
    }
    /* if datum has been set, check lost or not */
    if (sat) {
        slip=rtk->ssat[sat-1].slip[0]||rtk->ssat[sat-1].slip[1];
        for (i=0,j=-1;i<n;i++) {
            if (sat==obs[i].sat) { j=i; break; }
        }
        if (j<0||exc[j]||slip) {
            datum[sat-1]=0;
            datum[MAXSAT]=1;
        }
        return;
    }

    /* set sat with highest elevation as datum for the first epoch */
    for (i=0;i<n&&i<MAXOBS;i++) {
        if (exc[i]) continue;
        satel[i]=azel[i*2+1];
    }
    j=d_max(satel,n);
    sat=obs[j].sat;
    datum[sat-1]=1;
}
/* calculate antenna corrections and corrected measurements ------------------*/
extern void CalCorrAntObs(rtk_t *rtk, const obsd_t *obs, const nav_t *nav, 
    const double *rs, const double *rr, const double *azel, double *L, 
    double *P, double *Lc, double *Pc)
{
    const prcopt_t *opt=&rtk->opt;
    double dantr[NFREQ]={0},dants[NFREQ]={0};
    const int sat=obs->sat;

    /* satellite and receiver antenna model */
    satantpcv(rs,rr,nav->pcvs+sat-1,dants);
    antmodel(opt->pcvr,opt->antdel[0],azel,1,dantr);

    /* phase windup model */
    windupcorr(rtk->sol.time,rs,rr,&rtk->ssat[sat-1].phw);

    /* corrected phase and code measurements */
    corr_meas(obs,nav,azel,opt,dantr,dants,rtk->ssat[sat-1].phw,L,P,Lc,Pc);
}
/* update solution status ----------------------------------------------------*/
extern void userUdSolution(rtk_t *rtk, const obsd_t *obs, int n, int stat)
{
    const prcopt_t *opt=&rtk->opt;
    int i,j;

    /* test # of valid satellites */
    rtk->sol.ns=0;
    for (i=0;i<n&&i<MAXOBS;i++) {
        for (j=0;j<opt->nf;j++) {
            if (!rtk->ssat[obs[i].sat-1].vsat[j]) continue;
            rtk->ssat[obs[i].sat-1].lock[j]++;
            rtk->ssat[obs[i].sat-1].outc[j]=0;
            if (j==0) rtk->sol.ns++;
        }
    }
    rtk->sol.stat=rtk->sol.ns<MIN_NSAT_SOL?SOLQ_NONE:stat;

    if (rtk->sol.stat==SOLQ_FIX) {
        for (i=0;i<3;i++) {
            rtk->sol.rr[i]=rtk->xa[i];
            rtk->sol.qr[i]=(float)rtk->Pa[i+i*rtk->na];
        }
        rtk->sol.qr[3]=(float)rtk->Pa[1];
        rtk->sol.qr[4]=(float)rtk->Pa[1+2*rtk->na];
        rtk->sol.qr[5]=(float)rtk->Pa[2];
    }
    else {
        for (i=0;i<3;i++) {
            rtk->sol.rr[i]=rtk->x[i];
            rtk->sol.qr[i]=(float)rtk->P[i+i*rtk->nx];
        }
        rtk->sol.qr[3]=(float)rtk->P[1];
        rtk->sol.qr[4]=(float)rtk->P[2+rtk->nx];
        rtk->sol.qr[5]=(float)rtk->P[2];
        rtk->nfix=0;
    }
    //rtk->sol.dtr[0]=rtk->x[IC(0,opt)];
    //rtk->sol.dtr[1]=rtk->x[IC(1,opt)]-rtk->x[IC(0,opt)];

    for (i=0;i<n&&i<MAXOBS;i++) for (j=0;j<opt->nf;j++) {
        rtk->ssat[obs[i].sat-1].snr[j]=obs[i].SNR[j];
    }
    for (i=0;i<MAXSAT;i++) for (j=0;j<opt->nf;j++) {
        if (rtk->ssat[i].slip[j]&3) rtk->ssat[i].slipc[j]++;
        if (rtk->ssat[i].fix[j]==2&&stat!=SOLQ_FIX) rtk->ssat[i].fix[j]=1;
    }
}
/* initialize state and covariance -------------------------------------------*/
extern void userInitX_dp(rtk_t *rtk, double xi, double var, int i)
{
    int j;
    rtk->x[i]=xi;
    for (j=0;j<rtk->nx;j++) {
        rtk->P[i+j*rtk->nx]=rtk->P[j+i*rtk->nx]=i==j?var:0.0;
    }
}