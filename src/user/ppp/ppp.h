#ifndef _PPP_H
#define _PPP_H

#include "user/user.h"

/* number and index of states */
#if 0
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP(opt)     ((opt)->dynamics?9:3)
#define NC(opt)     (NSYS)
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt==TROPOPT_EST?1:3))
#define NI(opt)     ((opt)->ionoopt==IONOOPT_EST?MAXSAT:0)
#define ND(opt)     ((opt)->nf>=3?1:0)
#define NR(opt)     (NP(opt)+NC(opt)+NT(opt)+NI(opt)+ND(opt))
#define NB(opt)     (NF(opt)*MAXSAT)
#define NX(opt)     (NR(opt)+NB(opt))
#define IC(s,opt)   (NP(opt)+(s))
#define IT(opt)     (NP(opt)+NC(opt))
#define II(s,opt)   (NP(opt)+NC(opt)+NT(opt)+(s)-1)
#define ID(opt)     (NP(opt)+NC(opt)+NT(opt)+NI(opt))
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1)
#else
#define UCK(opt)    (((opt)->ionoopt==IONOOPT_EST&&(opt)->modear>ARMODE_OFF)?1:0)
#define PNF(opt)    ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define PNP(opt)    (3)
#define PNC(opt)    (NSYS)
#define PNT(opt)    ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt==TROPOPT_EST?1:3))
#define PNI(opt)    ((opt)->ionoopt==IONOOPT_EST?MAXSAT:0)
#define PND(opt)	(UCK(opt)?PNF(opt):0)
#define PNR(opt)    (PNP(opt)+PNC(opt)+PNT(opt)+PNI(opt)+PND(opt))
#define PNB(opt)    (PNF(opt)*MAXSAT)
#define PNX(opt)    (PNR(opt)+PNB(opt))

#define PIC(s,opt)   (PNP(opt)+(s-1)) /* s->sys (GPS:1) */
#define PIT(opt)     (PNP(opt)+PNC(opt))
#define PII(s,opt)	 (PNP(opt)+PNC(opt)+PNT(opt)+(s)-1)
#define PID(f,opt)	 (PNP(opt)+PNC(opt)+PNT(opt)+PNI(opt)+(f-1)) /* receiver phase bias */
#define PIB(s,f,opt) (PNR(opt)+(s-1)*PNF(opt)+(f-1))
#endif

#define TRACE_PPP	0

/* userPostpos_p.c -----------------------------------------------------------*/
extern int userPostPos_p(gtime_t ts, gtime_t te, double ti, double tu,
    const prcopt_t* popt, const solopt_t* sopt, const filopt_t* fopt);

/* userRTKPos_p.c ------------------------------------------------------------*/
extern int userRTKPos_p(rtk_t *rtk, obsd_t *obs, int n, const nav_t *nav);

/* userPPPos_p.c -------------------------------------------------------------*/
extern int userPPPnx_p(const prcopt_t *opt);
extern void userPPPos_p(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav);

/* userUdStates_p.c ----------------------------------------------------------*/
extern void userUdStates_p(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav);

/* userMWEst_cnes.c ----------------------------------------------------------*/
extern void userMWEst_cnes(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav, 
    const double *rs, const double *var_rs, const int *svh, const double *dr, 
    int *exc);

/* userMWAmbResol_cnes.c -----------------------------------------------------*/
extern int userMWAmbResol_cnes(rtk_t *rtk, const obsd_t *obs, int n, 
    const int *exc, double *x, double *P);
extern int userAmbResol_cnes2(rtk_t *rtk, const obsd_t *obs, int n, 
    const int *exc, const double *azel, double *x, double *P);

/* userPPPRes_p.c ------------------------------------------------------------*/
extern int userPPPRes_p(int post, const obsd_t *obs, int n, const double *rs,
    const double *dts, const double *var_rs, const double *var_dt, 
    const int *svh, const double *dr, int *exc, const nav_t *nav, 
    const double *xp, double *Pp, rtk_t *rtk, double *vl, double *HL, 
    double *DL, double *vp, double *HP, double *DP, double *azel, int *ix,
    int *nc);

/* userAmbResol_uup.c --------------------------------------------------------*/
extern int userAmbResol_cnes(rtk_t *rtk, const obsd_t *obs, int n, 
    const nav_t *nav, const int *exc, const double *azel, double *x, double *P);
extern int userAmbResol_uup(rtk_t *rtk, const obsd_t *obs, int n, 
    const nav_t *nav, const int *exc, const double *azel, double *x, double *P);

/* userOutSolStat_p.c --------------------------------------------------------*/
extern int userOpenSolStat_p(const char *file, int level);
extern void userCloseSolStat_p(void);
extern void userOutSolStat_p(rtk_t *rtk);

/* logRclk_p.c ---------------------------------------------------------------*/
extern int logRclkOpen_p(const char *file, const prcopt_t *opt);
extern void logRclkClose_p(void);
extern void logRclk_p(rtk_t *rtk);

/* logZtd_p.c ----------------------------------------------------------------*/
extern int logZtdOpen_p(const char *file, const prcopt_t *opt);
extern void logZtdClose_p(void);
extern void logZtd_p(rtk_t *rtk);

/* logIon_p.c ----------------------------------------------------------------*/
extern int logIonOpen_p(const char *file, const prcopt_t *opt);
extern void logIonClose_p(void);
extern void logIon_p(rtk_t *rtk, const obsd_t *obs, int n);

/* logAmb_p.c ----------------------------------------------------------------*/
extern int logAmbOpen_p(const char *file, const prcopt_t *opt);
extern void logAmbClose_p(void);
extern void logAmb_p(rtk_t *rtk, const obsd_t *obs, int n);

/* logRes_p.c ----------------------------------------------------------------*/
extern int logResOpen_p(const char *file, const prcopt_t *opt);
extern void logResClose_p(void);
extern void logRes_p(rtk_t *rtk, const obsd_t *obs, int n);

#endif /* _PPP_H */
