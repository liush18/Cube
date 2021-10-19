/* SBAS, ephemeris to satellite clock and position */
#include "sat.h"

/* sbas ephemeris to satellite clock bias --------------------------------------
* compute satellite clock bias with sbas ephemeris
* args  :gtime_t time     I   time by satellite clock (gpst)
*          seph_t *seph     I   sbas ephemeris
* return:satellite clock bias (s)
* notes :see ref [3]
*-----------------------------------------------------------------------------*/
extern double seph2clk(gtime_t time, const seph_t *seph)
{
    double t;
    int i;

    trace(4,"seph2clk: time=%s sat=%2d\n",time_str(time,3),seph->sat);

    t=timediff(time,seph->t0);

    for (i=0;i<2;i++) {
        t-=seph->af0+seph->af1*t;
    }
    return seph->af0+seph->af1*t;
}
/* sbas ephemeris to satellite position and clock bias -------------------------
* compute satellite position and clock bias with sbas ephemeris
* args  :gtime_t time     I   time (gpst)
*          seph_t  *seph    I   sbas ephemeris
*          double  *rs      O   satellite position {x,y,z} (ecef) (m)
*          double  *dts     O   satellite clock bias (s)
*          double  *var     O   satellite position and clock variance (m^2)
* return:none
* notes :see ref [3]
*-----------------------------------------------------------------------------*/
extern void seph2pos(gtime_t time, const seph_t *seph, double *rs, double *dts,
    double *var)
{
    double t;
    int i;

    trace(4,"seph2pos: time=%s sat=%2d\n",time_str(time,3),seph->sat);

    t=timediff(time,seph->t0);

    for (i=0;i<3;i++) {
        rs[i]=seph->pos[i]+seph->vel[i]*t+seph->acc[i]*t*t/2.0;
    }
    *dts=seph->af0+seph->af1*t;

    *var=var_uraeph(SYS_SBS,seph->sva);
}