#ifndef _SAT_H
#define _SAT_H

#include "base/base.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifdef WIN_DLL
#define EXPORT __declspec(dllexport) /* for Windows DLL */
#else
#define EXPORT
#endif /* WIN_DLL */

EXPORT void satantoff(gtime_t time, const double *rs, int sat, 
    const nav_t *nav, double *dant);
EXPORT int peph2pos(gtime_t time, int sat, const nav_t *nav, int opt,
    double *rs, double *dts, double *var);
EXPORT int peph_to_pos(gtime_t time, int sat, const nav_t *nav, int opt,
    double *rs, double *dts, double *vars, double *varc);

EXPORT double var_uraeph(int sys, int ura);
EXPORT double eph2clk(gtime_t time, const eph_t *eph);
EXPORT void eph2pos(gtime_t time, const eph_t *eph, double *rs, double *dts,
    double *var);
EXPORT double geph2clk(gtime_t time, const geph_t *geph);
EXPORT void geph2pos(gtime_t time, const geph_t *geph, double *rs, double *dts,
    double *var);
EXPORT double seph2clk(gtime_t time, const seph_t *seph);
EXPORT void seph2pos(gtime_t time, const seph_t *seph, double *rs, double *dts,
    double *var);
EXPORT int ephclk(gtime_t time, gtime_t teph, int sat, const nav_t *nav,
    double *dts);
EXPORT int ephpos(gtime_t time, gtime_t teph, int sat, const nav_t *nav,
    int iode, double *rs, double *dts, double *var, int *svh);
EXPORT int getseleph(int sys);

EXPORT void satposs(gtime_t teph, const obsd_t *obs, int n, const nav_t *nav,
    int ephopt, double *rs, double *dts, double *var, int *svh);
EXPORT void sat_poss(gtime_t teph, char *rcv, const obsd_t *obs, int n, 
    const nav_t *nav, double *rs, double *dts, double *vars, double *varc, 
    int *svh);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* _SAT_H */
