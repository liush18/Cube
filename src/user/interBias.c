
#include "user.h"

#define SQR(x)       ((x)*(x))
#define MAXDTE       900.0           /* max time difference to phase bias time (s) */
#define EXTERR_BIAS  1E-6            /* extrapolation error for bias (cycle/s) */

/* satellite carrier-phase bias interpolation ----------------------------------
* compute satellite carrier-phase bias and variance
* args   : prcopt_t  *opt      I   options
*          gtime_t   time      I   interpolation time
*          int       sat       I   satellite number
*          nav_t     *nav      I   navigation data
*          double    *bias     O   carrier-phase bias
*          double    *varBias  O   carrier-phase std
* return : (0:error, 1:ok)
* notes  : bias   [0:nf] = sat carrier-phase bias for nf frequencies (cycle)
*          varBias[0:nf] = sat carrier-phase bias std (cycle)
*-----------------------------------------------------------------------------*/
extern int interpol_phasebias(const prcopt_t *opt, gtime_t time, int sat, 
    const nav_t *nav, double *bias, double *varBias)
{
    double t[2],b[2],std;
    char str[32],id[8];
    int i,j,k,index;

    time2str(time,str,2);
    trace(4,"interpol_phasebias: time=%s sat=%2d\n",str,sat);

    if (nav->nb<2||
        timediff(time,nav->bias[0].time)<-MAXDTE||
        timediff(time,nav->bias[nav->nb-1].time)>MAXDTE) {
        satno2id(sat,id);
        trace(2,"no carrier-phase bias %s sat=%s\n",str,id);
        return 0;
    }
    /* binary search */
    for (i=0,j=nav->nb-1;i<j;) {
        k=(i+j)/2;
        if (timediff(nav->bias[k].time,time)<0.0) i=k+1; else j=k;
    }
    index=i<=0?0:i-1;

    /* linear interpolation for clock */
    t[0]=timediff(time,nav->bias[index  ].time);
    t[1]=timediff(time,nav->bias[index+1].time);

    /* carrier-phase bias for each frequency */
    for (i=0;i<opt->nf;i++) {

        bias[i]=0.0;
        varBias[i]=0.0;

        b[0]=nav->bias[index  ].bias[sat-1][i];
        b[1]=nav->bias[index+1].bias[sat-1][i];

        if (t[0]<=0.0) { /* before interpolation interval */
            if ((bias[i]=b[0])==0.0) return 0;
            std=nav->bias[index].std[sat-1][i]-EXTERR_BIAS*t[0];
        }
        else if (t[1]>=0.0) { /* after interpolation interval */
            if ((bias[i]=b[1])==0.0) return 0;
            std=nav->bias[index+1].std[sat-1][i]+EXTERR_BIAS*t[1];
        }
        else if (b[0]!=0.0&&b[1]!=0.0) {
            bias[i]=(b[1]*t[0]-b[0]*t[1])/(t[0]-t[1]);
            k=t[0]<-t[1]?0:1;   /* the nearest */
            std=nav->bias[index+k].std[sat-1][i]+EXTERR_BIAS*fabs(t[k]);
        }
        else {
            trace(2,"carrier-phase bias outage %s sat=%s\n",str,id);
            return 0;
        }
        varBias[i]=SQR(std);
    }
    return 1;
}

/* satellite code bias interpolation -------------------------------------------
* compute satellite code bias and variance
* args   : prcopt_t  *opt      I   options
*          gtime_t   time      I   interpolation time
*          int       sat       I   satellite number
*          nav_t     *nav      I   navigation data
*          double    *ucd      O   code bias drift
*          double    *varUcd   O   code bias drift std
* return : (0:error, 1:ok)
* notes  : ucd   [0:nf] = sat code bias for nf frequencies (cycle)
*          varUcd[0:nf] = sat code bias std (cycle)
*-----------------------------------------------------------------------------*/
extern int interpol_codebias(const prcopt_t *opt, gtime_t time, int sat, 
    const nav_t *nav, double *ucd, double *varUcd)
{
    double t[2],b[2],std;
    char str[32],id[8];
    int i,j,k,index;

    time2str(time,str,2);
    trace(4,"interpol_codebias: time=%s sat=%2d\n",str,sat);

    if (nav->nd<2||
        timediff(time,nav->ucd[0].time)<-MAXDTE||
        timediff(time,nav->ucd[nav->nd-1].time)>MAXDTE) {
        satno2id(sat,id);
        trace(2,"%s no ucd sat=%s\n",str,id);
        return -1;
    }
    /* binary search */
    for (i=0,j=nav->nd-1;i<j;) {
        k=(i+j)/2;
        if (timediff(nav->ucd[k].time,time)<0.0) i=k+1; else j=k;
    }
    index=i<=0?0:i-1;

    /* linear interpolation for clock */
    t[0]=timediff(time,nav->ucd[index  ].time);
    t[1]=timediff(time,nav->ucd[index+1].time);

    /* carrier-phase ucd for each frequency */
    for (i=0;i<opt->nf;i++) {

        ucd[i]=0.0;
        varUcd[i]=0.0;

        b[0]=nav->ucd[index  ].bias[sat-1][i];
        b[1]=nav->ucd[index+1].bias[sat-1][i];

        if (t[0]<=0.0) { /* before interpolation interval */
            if ((ucd[i]=b[0])==0.0) return 0;
            std=nav->ucd[index].std[sat-1][i]-EXTERR_BIAS*t[0];
        }
        else if (t[1]>=0.0) { /* after interpolation interval */
            if ((ucd[i]=b[1])==0.0) return 0;
            std=nav->ucd[index+1].std[sat-1][i]+EXTERR_BIAS*t[1];
        }
        else if (b[0]!=0.0&&b[1]!=0.0) {
            ucd[i]=(b[1]*t[0]-b[0]*t[1])/(t[0]-t[1]);
            k=t[0]<-t[1]?0:1;   /* the nearest */
            std=nav->ucd[index+k].std[sat-1][i]+EXTERR_BIAS*fabs(t[k]);
        }
        else {
            trace(2,"ucd outage %s sat=%s\n",str,id);
            return -1;
        }
        varUcd[i]=SQR(std);
    }
    return 1;
}
