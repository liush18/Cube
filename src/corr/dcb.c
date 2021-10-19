
#include "corr.h"

#define SQR(x)      ((x)*(x))


/* antenna corrected measurements --------------------------------------------*/
extern int corrObsMeasurement(const obsd_t *obs, const nav_t *nav, const double *azel,
    const prcopt_t *opt, const double *dantr, const double *dants, double phw, 
    double *L, double *P, double *Lc, double *Pc)
{
    double freq[NFREQ]={0},C1,C2,P1_P2=30.0;
    int i,sys=satsys(obs->sat,NULL);

    for (i=0;i<NFREQ;i++) {
        L[i]=P[i]=0.0;
        freq[i]=sat2freq(obs->sat,obs->code[i],nav);
        if (freq[i]==0.0||obs->L[i]==0.0||obs->P[i]==0.0) continue;
        if (testsnr(0,0,azel[1],obs->SNR[i]*SNR_UNIT,&opt->snrmask)) continue;

        /* antenna phase center and phase windup correction */
        L[i]=obs->L[i]*CLIGHT/freq[i]-dants[i]-dantr[i]-phw*CLIGHT/freq[i];
        P[i]=obs->P[i]-dants[i]-dantr[i];

        /* P1-C1,P2-C2 dcb correction (C1->P1,C2->P2) */
        if (sys==SYS_GPS||sys==SYS_GLO) {
            if (obs->code[i]==CODE_L1C) P[i]+=nav->cbias[obs->sat-1][1];
            if (obs->code[i]==CODE_L2C) P[i]+=nav->cbias[obs->sat-1][2];
        }
    }
    if (fabs(P[0]-P[1])>P1_P2) {
        trace(2,"%s reject sat=%d by P1_P2 detection\n",
            time_str(obs->time,2),obs->sat);
        return 0;
    }
    /* iono-free LC */
    if (freq[0]==0.0||freq[1]==0.0) return 0;
    C1= SQR(freq[0])/(SQR(freq[0])-SQR(freq[1]));
    C2=-SQR(freq[1])/(SQR(freq[0])-SQR(freq[1]));
    if (L[0]!=0.0&&L[1]!=0.0) *Lc=C1*L[0]+C2*L[1];
    if (P[0]!=0.0&&P[1]!=0.0) *Pc=C1*P[0]+C2*P[1];

    return 1;
}
/* antenna corrected measurements --------------------------------------------*/
extern void corr_meas(const obsd_t *obs, const nav_t *nav, const double *azel,
    const prcopt_t *opt, const double *dantr, const double *dants, double phw,
    double *L, double *P, double *Lc, double *Pc)
{
    double freq[NFREQ]={0},C1,C2;
    int i,sys=satsys(obs->sat,NULL);

    for (i=0;i<NFREQ;i++) {
        L[i]=P[i]=0.0;
        freq[i]=sat2freq(obs->sat,obs->code[i],nav);
        if (freq[i]==0.0||obs->L[i]==0.0||obs->P[i]==0.0) continue;
        if (testsnr(0,0,azel[1],obs->SNR[i]*SNR_UNIT,&opt->snrmask)) continue;

        /* antenna phase center and phase windup correction */
        L[i]=obs->L[i]*CLIGHT/freq[i]-dants[i]-dantr[i]-phw*CLIGHT/freq[i];
        P[i]=obs->P[i]-dants[i]-dantr[i];

        /* P1-C1,P2-C2 dcb correction (C1->P1,C2->P2) */
        if (sys==SYS_GPS||sys==SYS_GLO) {
            if (obs->code[i]==CODE_L1C) P[i]+=nav->cbias[obs->sat-1][1];
            if (obs->code[i]==CODE_L2C) P[i]+=nav->cbias[obs->sat-1][2];
        }
    }
    /* iono-free LC */
    *Lc=*Pc=0.0;
    if (freq[0]==0.0||freq[1]==0.0) return;
    C1= SQR(freq[0])/(SQR(freq[0])-SQR(freq[1]));
    C2=-SQR(freq[1])/(SQR(freq[0])-SQR(freq[1]));

    if (L[0]!=0.0&&L[1]!=0.0) *Lc=C1*L[0]+C2*L[1];
    if (P[0]!=0.0&&P[1]!=0.0) *Pc=C1*P[0]+C2*P[1];
}