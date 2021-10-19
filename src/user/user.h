#ifndef _USER_H
#define _USER_H

#include "corr/corr.h"
#include "io/io.h"
#include "qc/qc.h"
#include "sat/sat.h"

typedef struct {        /* solution type */
    gtime_t time;       /* time (GPST) */
    double rr[6];       /* position/velocity (m|m/s) */
                        /* {x,y,z,vx,vy,vz} or {e,n,u,ve,vn,vu} */
    float  qr[6];       /* position variance/covariance (m^2) */
                        /* {c_xx,c_yy,c_zz,c_xy,c_yz,c_zx} or */
                        /* {c_ee,c_nn,c_uu,c_en,c_nu,c_ue} */
    float  qv[6];       /* velocity variance/covariance (m^2/s^2) */
    double dtr[6];      /* receiver clock bias to time systems (s) */
    unsigned char type; /* type (0:xyz-ecef,1:enu-baseline) */
    unsigned char stat; /* solution status (SOLQ_???) */
    unsigned char ns;   /* number of valid satellites */
    float age;          /* age of differential (s) */
    float ratio;        /* AR ratio factor for valiation */
    float rate;         /* AR success rate */
    float thres;        /* AR ratio threshold for valiation */
} sol_t;

typedef struct {        /* ambiguity control type */
    gtime_t epoch[4];   /* last epoch */
    int n[4];           /* number of epochs */
    double LC [4];      /* linear combination average */
    double LCv[4];      /* linear combination variance */
    int fixcnt;         /* fix count */
    char flags[MAXSAT]; /* fix flags */
} ambc_t;

typedef struct {
    double floamb[NFREQ];   /* float ambiguity, sd for if {wide,narr} */
    double flostd[NFREQ];
    int fixamb[NFREQ];      /* fixed ambiguity, sd for if {wide,narr} */
    int fixflag;
}ambr_t;

typedef struct {        /* RTK control/result type */
    sol_t  sol;         /* RTK solution */
    double rb[6];       /* base position/velocity (ecef) (m|m/s) */
    int nx,na,nw;       /* number of float states/fixed states */
    double tt;          /* time difference between current and previous (s) */
    double *x,*P;       /* float states and their covariance */
    double *xa,*Pa;     /* fixed states and their covariance */
    double *xw,*Pw;     /* wide-lane ambiguity estimation */
    int nfix;           /* number of continuous fixes of ambiguity */
    ambc_t ambc[MAXSAT];/* ambibuity control */
    ssat_t ssat[MAXSAT];/* satellite status */
    int neb;            /* bytes in error message buffer */
    char errbuf[MAXERRMSG]; /* error message buffer */
    prcopt_t opt;       /* processing options */
    int jumpc;          /* clock jump values --by ls */
    int datum[MAXSAT+1];/* decoupled model datum of ambiguities */
                        /* last index means datum has been set */
    ambt_t amb[MAXSAT]; /* ar tracking */
    ambr_t ambr[MAXSAT];/* record ambiguity */
    char cTime[32];     /* time string of current epoch */
} rtk_t;

/* pntpos.c ------------------------------------------------------------------*/
extern int ionocorr(gtime_t time, const nav_t *nav, int sat, const double *pos,
    const double *azel, int ionoopt, double *ion, double *var);
extern int tropcorr(gtime_t time, const nav_t *nav, const double *pos,
    const double *azel, int tropopt, double *trp, double *var);
extern int pntpos(const obsd_t *obs, int n, const nav_t *nav,
    const prcopt_t *opt, sol_t *sol, double *azel, ssat_t *ssat,
    char *msg);

/* interBias.c ---------------------------------------------------------------*/
extern int interpol_phasebias(const prcopt_t *opt, gtime_t time, int sat, 
    const nav_t *nav, double *bias, double *varBias);
extern int interpol_codebias(const prcopt_t *opt, gtime_t time, int sat, 
    const nav_t *nav, double *ucd, double *varUcd);

/* solution.c ----------------------------------------------------------------*/
extern int outprcopts(unsigned char *buff, const prcopt_t *opt);
extern int outsolheads(unsigned char *buff, const solopt_t *opt);
extern int outsols(unsigned char *buff, const sol_t *sol, const double *rb,
    const solopt_t *opt);
extern void outprcopt(FILE *fp, const prcopt_t *opt);
extern void outsolhead(FILE *fp, const solopt_t *opt);
extern void outsol(FILE *fp, const sol_t *sol, const double *rb,
    const solopt_t *opt);

/* userSession.c -------------------------------------------------------------*/
extern int userOpenSession(const prcopt_t *popt, const filopt_t *fopt, 
    nav_t *nav, pcvs_t *pcv);
extern void userCloseSession(nav_t *nav,pcvs_t *pcv);

/* userCmn.c -----------------------------------------------------------------*/
extern pcvs_t *userGetPcv();
extern obs_t  *userGetObs();
extern nav_t  *userGetNav();
extern sta_t  *userGetSta();
extern int userGet_revs();
extern int userGet_aborts();
extern int userJudge_sol();
extern void userSet_revs(int i);
extern void userSet_iobs(int i);
extern void userSet_isol(int i);
extern void userResetPar();
extern void userFree_sol_rb();
extern void userMallocSol();
extern int userCombForw(rtk_t *rtk);
extern int userCombBack(rtk_t *rtk);
extern void userSortObs();
extern int userCheckbrk(const char *format,...);
extern int userInputObs(obsd_t *obs, int solq, const prcopt_t *popt);
extern void userCombres(FILE *fp,const prcopt_t *popt,const solopt_t *sopt);
extern void userFreeObs();
extern int UserOutHead(const char *outfile, const prcopt_t *popt,
    const solopt_t *sopt);
extern void userRTKInit(rtk_t *rtk, const prcopt_t *opt, int nx);
extern void userRTKFree(rtk_t *rtk);

/* userCmnPPP.c --------------------------------------------------------------*/
extern void userSetAmbDatum(const rtk_t *rtk, const obsd_t *obs, int n, 
    const int *exc, const double *azel, int *datum);
extern void CalCorrAntObs(rtk_t *rtk, const obsd_t *obs, const nav_t *nav, 
    const double *rs, const double *rr, const double *azel, double *L, 
    double *P, double *Lc, double *Pc);
extern void userUdSolution(rtk_t *rtk, const obsd_t *obs, int n, int stat);
extern void userInitX_dp(rtk_t *rtk, double xi, double var, int i);

#endif /* _USER_H */
