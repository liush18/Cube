#ifndef _CORR_H
#define _CORR_h

#include "base\base.h"

/* dcb.c */
extern int corrObsMeasurement(const obsd_t *obs, const nav_t *nav, const double *azel,
    const prcopt_t *opt, const double *dantr, const double *dants, double phw, 
    double *L, double *P, double *Lc, double *Pc);
extern void corr_meas(const obsd_t *obs, const nav_t *nav, const double *azel,
    const prcopt_t *opt, const double *dantr, const double *dants, double phw,
    double *L, double *P, double *Lc, double *Pc);

/* pcv.c */
extern void antmodel(const pcv_t* pcv, const double* del, const double* azel,
    int opt, double* dant);
extern void satantpcv(const double* rs, const double* rr, const pcv_t* pcv,
    double* dant);
extern pcv_t* searchpcv(int sat, const char *type, gtime_t time,
    const pcvs_t *pcvs);
extern void posSetpcv(gtime_t time, prcopt_t *popt, nav_t *nav, 
    const pcvs_t *pcvs, const pcvs_t *pcvr, const sta_t *sta);

/* phw.c */
extern int model_phw(gtime_t time, int sat, const char* type, int opt,
    const double* rs, const double* rr, double* phw);
extern void windupcorr(gtime_t time, const double *rs, const double *rr,
    double *phw);

/* trop.c */
extern double tropmodel(gtime_t time, const double *pos, const double* azel,
    double humi);
extern double tropMops(gtime_t time, const double *pos, const double *azel,
    double *var);
extern double tropNMF(gtime_t time, const double *pos, const double *azel,
    double *mapfw);
extern void tropModel(gtime_t time, const double *pos, const double *azel,
    const double *x, double *dtdx, double *dtrp, double *var);
//extern double tropModel2(gtime_t time, const double *pos, const double *azel,
//    double *m_w, double *zhd, double *var);

/* iono.c */
extern double ionmodel(gtime_t t, const double *ion, const double *pos,
    const double *azel);
extern double ionmapf(const double *pos, const double *azel);
extern double ionppp(const double *pos, const double *azel, double re,
    double hion, double *posp);
extern int iontec(gtime_t time, const nav_t *nav, const double *pos,
    const double *azel, int opt, double *delay, double *var);

/* tide.c */
extern int geterp(const erp_t* erp, gtime_t time, double* erpv);
extern void tidedisp(gtime_t tutc, const double *rr, int opt, const erp_t *erp,
    const double *odisp, double *dr);


#endif /* _CORR_H */
