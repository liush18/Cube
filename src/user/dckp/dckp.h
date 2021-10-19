#ifndef _DCKP_H
#define _DCKP_H

#include "user/user.h"

/* number and index of states */
#define DPNF(opt)	  (2)	    /* ambiguity number of dck model */
#define DPCK(opt)	  (3)       /* clock number of dck model */
#define DPNP(opt)     (3)       /* position */
#define DPNC(opt)     (DPCK(opt)*NSYS)                  /* clocks */
#define DPNT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt==TROPOPT_EST?1:3)) /* troposphere */
#define DPNR(opt)     (DPNP(opt)+DPNC(opt)+DPNT(opt))   /* number for receiver */
#define DPNB(opt)     (DPNF(opt)*MAXSAT)                /* ambiguity */
#define DPNX(opt)     (DPNR(opt)+DPNB(opt))             /* total */

#define DPIC(sys,n,opt)     (DPNP(opt)+(sys-1)*DPCK(opt)+(n-1)) /* sys (GPS:1) */
#define DPIT(opt)           (DPNP(opt)+DPNC(opt))
#define DPIB(s,f,opt)       (DPNR(opt)+(s-1)*DPNF(opt)+(f-1))   /* 1:narr, 2:wide */


#define TRACE_PPP	0


/* userPostPos_dp.c ----------------------------------------------------------*/
extern int userPostPos_dp(gtime_t ts, gtime_t te, double ti, double tu,
    const prcopt_t *popt, const solopt_t *sopt, const filopt_t *fopt);

/* userRTKPos_dp -------------------------------------------------------------*/
extern int userRTKPos_dp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav);

/* userPPPos_dp.c ------------------------------------------------------------*/
extern int userPPPnx_dp(const prcopt_t *opt);
extern void userPPPos_dp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav);

/* userPPPRes_dp -------------------------------------------------------------*/
extern int userPPPRes_dp(int post, const obsd_t *obs, int n, const double *rs,
    const double *dts, const double *var_rs, const double *var_dt, 
    const int *svh, const double *dr, int *exc, const nav_t *nav, 
    const double *xp, double *Pp, rtk_t *rtk, double *vl, double *HL, 
    double *DL, double *vp, double *HP, double *DP, double *azel, int *ix, 
    int *nc);

/* userUdStates_dp.c ---------------------------------------------------------*/
extern void userUdStates_dp(rtk_t *rtk, const obsd_t *obs, int n, 
    const nav_t *nav);

/* userAmbResol_dp.c ---------------------------------------------------------*/
extern int userAmbResol_dp(rtk_t *rtk, const obsd_t *obs, int n, const int *exc,
    const double *azel, double *x, double *P);

/* logRclk_dp.c --------------------------------------------------------------*/
extern int logRclkOpen_dp(const char *file, const prcopt_t *opt);
extern void logRclkClose_dp(void);
extern void logRclk_dp(rtk_t *rtk);

/* logZtd_dp.c ---------------------------------------------------------------*/
extern int logZtdOpen_dp(const char *file, const prcopt_t *opt);
extern void logZtdClose_dp(void);
extern void logZtd_dp(rtk_t *rtk);

/* logAmb_dp.c ---------------------------------------------------------------*/
extern int logAmbOpen_dp(const char *file, const prcopt_t *opt);
extern void logAmbClose_dp(void);
extern void logAmb_dp(rtk_t *rtk, const obsd_t *obs, int n);

/* logRes_dp.c ---------------------------------------------------------------*/
extern int logResOpen_dp(const char *file, const prcopt_t *opt);
extern void logResClose_dp(void);
extern void logRes_dp(rtk_t *rtk, const obsd_t *obs, int n);


#endif /* _DCKP_H */
