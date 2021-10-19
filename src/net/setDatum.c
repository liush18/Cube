
#include "net.h"

#define DTMELE          20.0*D2R    /* sat elevation mask of datum */
#define SQR(x)		    (x*x)

#define EPOCH_RELOCK1    0
#define EPOCH_RELOCK2    120         /* re-lock satellite */

/* set SSC to rcv with highest elevation -------------------------------------*/
static void setSSC(const char *str, rcv_t *rcv, int n, int sat, int *SSC,
    const int *SRC)
{
    double satel[MAXRCV];
    int i;

    /* current sat elevation in each rcv */
    for (i=0;i<n&&i<MAXRCV;i++) {
        /* if sat used to estimate rcv, not set as SSC (after init) */
        if (rcv[i].nobs==0||rcv[i].datum[sat-1]<OBS_USE||SRC[i]==sat) {
            satel[i]=0.0;
        }
        else satel[i]=rcv[i].ssat[sat-1].azel[1];
        if (satel[i]<DTMELE) satel[i]=0.0;
    }
    /* rcv with highest elevation of current sat */
    if ((i=d_max(satel,n))>=0) {
        SSC[sat-1]=i+1;
        trace(2,"%s add new sat(%d) as SSC to rcv(%s)\n",str,sat,rcv[i].sta.name);
    }
    else {
        SSC[sat-1]=0; netRejSatAllRcv(rcv,n,sat);
        trace(2,"%s could not add sat(%d) as SSC to rcv\n",str,sat);
    }
}
/* set SRC for rcv ircv ------------------------------------------------------*/
static int setRcvSRC(const rcv_t *rcv, int ircv, const int *SSC, int *SRC)
{
    const ssat_t *ssat=rcv[ircv-1].ssat;
    int i,j,sat;

    /* measurement observed by rcv ircv */
    for (i=j=0;i<rcv[ircv-1].nobs&&i<MAXOBS;i++) {
        sat=rcv[ircv-1].obs[i].sat;
        if (!SSC[sat-1]) continue;  /* not estimable satellite clock */
        if (rcv[ircv-1].datum[sat-1]==OBS_REJ) continue;
        if (ssat[sat-1].azel[1]<DTMELE) continue;
        if (!j||ssat[sat-1].azel[1]>ssat[j-1].azel[1]) j=sat;
    }
    /* satellite with highest elevation as SRC */
    if (j) { SRC[ircv-1]=j; return 1; }

    return 0;
}
/* set SSC of rcv ircv -------------------------------------------------------*/
static void setRcvSSC(const rcv_t *rcv, int ircv, int *SSC)
{
    int i,sat;
    for (i=0;i<rcv[ircv-1].nobs&&i<MAXOBS;i++) {
        sat=rcv[ircv-1].obs[i].sat;
        if (SSC[sat-1]) continue;
        if (rcv[ircv-1].datum[sat-1]==OBS_REJ) continue;
        if (rcv[ircv-1].ssat[sat-1].azel[1]<DTMELE) continue;
        SSC[sat-1]=ircv;  /* clock of sat been separated */
    }
}
/* initializing ambiguity datum ----------------------------------------------*/
static void initDatum(const net_t *uck, rcv_t *rcv, int n, int *SSC, int *SRC)
{
    const int ref=uck->opt.ircv;
    int i;

    /* SSC of reference rcv */
    setRcvSSC(rcv,ref,SSC);

    /* SRC determining */
    for (i=0;i<n&&i<MAXRCV;i++) {
        if (i+1==ref) continue;
        /* processing rcv i+1 */
        if (setRcvSRC(rcv,i+1,SSC,SRC)) {   /* set SRC */
            setRcvSSC(rcv,i+1,SSC);         /* set SSC */
        }
        else {
            trace(2,"%s could not set SRC for rcv(%s) nobs=%d\n",
                uck->time,rcv[i].sta.name,rcv[i].nobs);
            rcv[i].nobs=0;
        }
    }
    /* processing satellite clock not estimable */
    for (i=0;i<MAXSAT;i++) {
        if (!SSC[i]) setSSC(uck->time,rcv,n,i+1,SSC,SRC);
    }
}
/* check datum ---------------------------------------------------------------*/
static void checkDatum(const net_t *net, rcv_t *rcv, int n, int *SSC, int *SRC)
{
    const int nx=net->nx;
    char id[4];
    int i,j,k,sat,ircv,slip;

    /* lost SRC */
    for (i=0;i<n&&i<MAXRCV;i++) {
        /* invalid rcv or out rcv in previous epoch */
        /* all obs objected -> maybe nobs!=0 but out in previous epoch */
        if (rcv[i].nobs==0||rcv[i].outc) { 
            /* remove SRC of invalid rcv */
            if (SRC[i]) {
                satno2id(SRC[i],id); SRC[i]=0;
                trace(3,"%s remove SRC(%s) of invalid rcv(%s)\n",
                    net->time,id,rcv[i].sta.name);
            }
            /* remove SSC of invalid rcv */
            for (j=0;j<MAXSAT;j++) {
                if (SSC[j]!=i+1) continue;
                satno2id(j+1,id); SSC[j]=0;
                trace(3,"%s remove SSC(%s) of invalid rcv(%s)\n",
                    net->time,id,rcv[i].sta.name);
            }
            continue; 
        }
        /* SRC of current rcv has been removed */
        if (!(sat=SRC[i])) continue;
        slip=rcv[i].ssat[sat-1].slip[0]||rcv[i].ssat[sat-1].slip[1];
        /* if slip or invalid obs */
        if (slip||rcv[i].datum[sat-1]<OBS_USE) {
            SRC[i]=0;
            satno2id(sat,id);
            trace(3,"%s remove SRC(%s) of rcv(%s) slip=%d outc=%d datum=%d\n",
                net->time,id,rcv[i].sta.name,slip,rcv[i].ssat[sat-1].outc[0],
                rcv[i].datum[sat-1]);
        }
    }

    /* lost SSC */
    for (i=0;i<MAXSAT;i++) {
        /* out sat in previous epoch */
        if (net->outc[i]) {
            /* remove invalid sat as SSC */
            if (SSC[i]) {
                satno2id(i+1,id); SSC[i]=0;
                trace(3,"%s remove invalid sat(%s) as SSC to rcv(%s)\n",
                    net->time,id,rcv[SSC[i]-1].sta.name,id);
            }
            /* remove invalid sat as SRC to rcv */
            for (j=0;j<n&&j<MAXRCV;j++) {
                if (SRC[j]!=i+1) continue;
                satno2id(i+1,id); SRC[j]=0;
                trace(3,"%s remove invalid sat(%s) as SRC to rcv(%s)\n",
                    net->time,id,rcv[j].sta.name,id);
            }
        }
        if (!(ircv=SSC[i])) continue;
        slip=rcv[ircv-1].ssat[i].slip[0]||rcv[ircv-1].ssat[i].slip[1];
        /* if slip or invalid obs */
        if (slip||rcv[ircv-1].datum[i]<OBS_USE) {
            SSC[i]=0;
            satno2id(i+1,id);
            trace(3,"%s remove SSC(%s) of rcv(%s) slip=%d outc=%d datum=%d\n",
                net->time,id,rcv[ircv-1].sta.name,slip,rcv[ircv-1].ssat[i].outc[0],
                rcv[ircv-1].datum[i]);
        }
    }

    /* new rcv */
    for (i=0;i<n&&i<MAXRCV;i++) {
        if (!rcv[i].outc||rcv[i].nobs==0) continue;
        if (rcv[i].outc>EPOCH_RELOCK1&&
            (rcv[i].outc%EPOCH_RELOCK2)>EPOCH_RELOCK1) {
            netRejRcvAllSat(&rcv[i]);
            continue;
        }
        /* for rcv unlock in previous epoch and lock in this epoch */
        for (j=k=0;j<rcv[i].nobs&&j<MAXOBS;j++) {
            sat=rcv[i].obs[j].sat;
            if (rcv[i].datum[sat-1]==OBS_REJ) continue;
            if (rcv[i].ssat[sat-1].azel[1]<DTMELE) continue;
            if (!k||rcv[i].ssat[sat-1].azel[1]>rcv[i].ssat[k-1].azel[1]) k=sat;
        }
        if (k) {
            SRC[i]=k;
            satno2id(k,id);
            trace(2,"%s add new rcv(%s) SRC(%s)\n",net->time,rcv[i].sta.name,id);
        }
        else {
            netRejRcvAllSat(&rcv[i]);
            trace(2,"%s could not add new rcv(%s)\n",net->time,rcv[i].sta.name);
        }
    }

    /* add new sat 
    if sat unlock by all rcv in previous epoch, try to add this epoch */
    for (i=0;i<MAXSAT;i++) {
        if (!net->outc[i]) continue;
        /* reject sat if outlock special times */
        if (net->outc[i]>EPOCH_RELOCK1&&
            (net->outc[i]%EPOCH_RELOCK2)>EPOCH_RELOCK1) {
            netRejSatAllRcv(rcv,n,i+1);
            continue;
        }
        /* need to set datum for this sat */
        setSSC(net->time,rcv,n,i+1,SSC,SRC);
    }
}
/* set ambiguity and update ambiguity ----------------------------------------*/
extern void setDatum(const net_t *uck, rcv_t *rcv, int n, int flag, int *SSC, 
    int *SRC)
{
    int i;

    /* ambiguity datum has not been initialized */
    if (!flag) initDatum(uck,rcv,n,SSC,SRC);
    else checkDatum(uck,rcv,n,SSC,SRC);

    /* update ambiguity datum */
    for (i=0;i<n&&i<MAXRCV;i++) {
        if (SRC[i]) rcv[i].datum[SRC[i]-1]=OBS_SRC;
    }
    for (i=0;i<MAXSAT;i++) {
        if (SSC[i]) rcv[SSC[i]-1].datum[i]=OBS_SSC;
    }
}