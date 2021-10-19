#ifndef _IO_H
#define _IO_H

#include "base/base.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifdef WIN_DLL
#define EXPORT __declspec(dllexport) /* for Windows DLL */
#else
#define EXPORT
#endif /* WIN_DLL */

#define NAVEXP      "D"                 /* exponent letter in RINEX NAV */
#define NUMSYS      7                   /* number of systems */
#define MAXRNXLEN   (16*MAXOBSTYPE+4)   /* max RINEX record length */
#define MAXPOSHEAD  1024                /* max head line position */
#define MINFREQ_GLO -7                  /* min frequency number GLONASS */
#define MAXFREQ_GLO 13                  /* max frequency number GLONASS */
#define NINCOBS     262144              /* incremental number of obs data */


/* options.c */
EXPORT opt_t sysopts[];              /* system options table */
EXPORT opt_t *searchopt(const char *name, const opt_t *opts);
EXPORT int str2opt(opt_t *opt, const char *str);
EXPORT int opt2str(const opt_t *opt, char *str);
EXPORT int opt2buf(const opt_t *opt, char *buff);
EXPORT int loadopts(const char *file, opt_t *opts);
EXPORT int saveopts(const char *file, const char *mode, const char *comment,
    const opt_t *opts);
EXPORT void resetsysopts(void);
EXPORT void getsysopts(prcopt_t *popt, solopt_t *sopt, filopt_t *fopt);
EXPORT void setsysopts(const prcopt_t *prcopt, const solopt_t *solopt,
    const filopt_t *filopt);


/* readblq.c */
EXPORT int readblq(const char* file, const char* sta, double* odisp);
EXPORT void readotl_pos(prcopt_t *popt, const char *file, const sta_t *sta);

/* readclk.c */
//EXPORT int readrnxc(const char* file, const int mask, nav_t* nav);

/* readdcb.c */
EXPORT int readdcb(const char* file, nav_t* nav, const sta_t* sta);

/* readerp.c */
EXPORT int readerp(const char* file, erp_t* erp);

/* readpcv.c */
EXPORT int readpcv(const char* file, pcvs_t* pcvs);

/* readsnx.c */


/* readsp3.c -----------------------------------------------------------------*/
EXPORT void combpeph(nav_t* nav, int opt);
EXPORT void readsp3(const char *file, nav_t *nav, int opt);

/* readtec.c -----------------------------------------------------------------*/
EXPORT int readtec(const char *file, nav_t *nav, int opt);


/* readrnx.c -----------------------------------------------------------------*/
EXPORT void uniqnav(nav_t *nav);
EXPORT int sortobs(obs_t *obs);
EXPORT int screent(gtime_t time, gtime_t ts, gtime_t te, double tint);
EXPORT int set_sysmask(const char *opt);
EXPORT int readrnxfile(const char *file, gtime_t ts, gtime_t te, double tint,
    const char *opt, int flag, int index, char *type, obs_t *obs, nav_t *nav, 
    sta_t *sta);
EXPORT int readrnxt(const char *file, int rcv, gtime_t ts, gtime_t te,
    double tint, const char *opt, obs_t *obs, nav_t *nav, sta_t *sta);
EXPORT int readrnx(const char *file, int rcv, const char *opt, obs_t *obs,
    nav_t *nav, sta_t *sta);
EXPORT int readrnxh(FILE *fp, double *ver, char *type, int *sys, int *tsys,
    char tobs[][MAXOBSTYPE][4], nav_t *nav, sta_t *sta);

/* readrnx_nav.c -------------------------------------------------------------*/
EXPORT void decode_navh(char *buff, nav_t *nav);
EXPORT void decode_gnavh(char *buff, nav_t *nav);
EXPORT void decode_hnavh(char *buff, nav_t *nav);
EXPORT int readrnxnav(FILE *fp, const char *opt, double ver, int sys, nav_t *nav);

/* readrnx_obs.c -------------------------------------------------------------*/
EXPORT void decode_obsh(FILE *fp, char *buff, double ver, int *tsys,
    char tobs[][MAXOBSTYPE][4], nav_t *nav, sta_t *sta);
EXPORT int readrnxobs(FILE *fp, gtime_t ts, gtime_t te, double tint,
    const char *opt, int rcv, double ver, int *tsys, char tobs[][MAXOBSTYPE][4],
    obs_t *obs, sta_t *sta);
EXPORT int readrnxobsb(FILE *fp, const char *opt, double ver, int *tsys,
    char tobs[][MAXOBSTYPE][4], int *flag, obsd_t *data, sta_t *sta);
EXPORT void saveslips(uint8_t slips[][NFREQ+NEXOBS], obsd_t *data);
EXPORT void restslips(uint8_t slips[][NFREQ+NEXOBS], obsd_t *data);

/* readrnx_clk.c -------------------------------------------------------------*/
EXPORT int readrnxclk(FILE *fp, const char *opt, int index, nav_t *nav);
EXPORT int readrnxc(const char *file, nav_t *nav);

/* readdck.c -----------------------------------------------------------------*/
EXPORT int readdck(const char *file, const int mask, nav_t *nav);

/* readBias.c ----------------------------------------------------------------*/
EXPORT int readBias(const char *file, const int mask, nav_t *nav);

/* readcob.c -----------------------------------------------------------------*/
EXPORT int readcob(const char *file, int mask, int nf, nav_t *nav);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* _IO_H */
