/* from ephemeris to space and time parameters of satellites */
#include "sat.h"

/* ephemeris selections ------------------------------------------------------*/
static int eph_sel[]={ /* GPS,GLO,GAL,QZS,BDS,IRN,SBS */
    0,0,0,0,0,0,0
};

/* select ephememeris --------------------------------------------------------*/
static eph_t *seleph(gtime_t time, int sat, int iode, const nav_t *nav)
{
    double t,tmax,tmin;
    int i,j=-1,sys,sel;

    trace(4,"seleph  : time=%s sat=%2d iode=%d\n",time_str(time,3),sat,iode);

    sys=satsys(sat,NULL);
    switch (sys) {
    case SYS_GPS: tmax=MAXDTOE+1.0    ; sel=eph_sel[0]; break;
    case SYS_GAL: tmax=MAXDTOE_GAL    ; sel=eph_sel[2]; break;
    case SYS_QZS: tmax=MAXDTOE_QZS+1.0; sel=eph_sel[3]; break;
    case SYS_CMP: tmax=MAXDTOE_CMP+1.0; sel=eph_sel[4]; break;
    case SYS_IRN: tmax=MAXDTOE_IRN+1.0; sel=eph_sel[5]; break;
    default: tmax=MAXDTOE+1.0; break;
    }
    tmin=tmax+1.0;

    for (i=0;i<nav->n;i++) {
        if (nav->eph[i].sat!=sat) continue;
        if (iode>=0&&nav->eph[i].iode!=iode) continue;
        if (sys==SYS_GAL) {
            sel=getseleph(SYS_GAL);
            if (sel==0&&!(nav->eph[i].code&(1<<9))) continue; /* I/NAV */
            if (sel==1&&!(nav->eph[i].code&(1<<8))) continue; /* F/NAV */
            if (timediff(nav->eph[i].toe,time)>=0.0) continue; /* AOD<=0 */
        }
        if ((t=fabs(timediff(nav->eph[i].toe,time)))>tmax) continue;
        if (iode>=0) return nav->eph+i;
        if (t<=tmin) {j=i; tmin=t;} /* toe closest to time */
    }
    if (iode>=0||j<0) {
        trace(3,"no broadcast ephemeris: %s sat=%2d iode=%3d\n",
            time_str(time,0),sat,iode);
        return NULL;
    }
    return nav->eph+j;
}
/* select glonass ephememeris ------------------------------------------------*/
static geph_t *selgeph(gtime_t time, int sat, int iode, const nav_t *nav)
{
    double t,tmax=MAXDTOE_GLO,tmin=tmax+1.0;
    int i,j=-1;

    trace(4,"selgeph : time=%s sat=%2d iode=%2d\n",time_str(time,3),sat,iode);

    for (i=0;i<nav->ng;i++) {
        if (nav->geph[i].sat!=sat) continue;
        if (iode>=0&&nav->geph[i].iode!=iode) continue;
        if ((t=fabs(timediff(nav->geph[i].toe,time)))>tmax) continue;
        if (iode>=0) return nav->geph+i;
        if (t<=tmin) {j=i; tmin=t;} /* toe closest to time */
    }
    if (iode>=0||j<0) {
        trace(3,"no glonass ephemeris  : %s sat=%2d iode=%2d\n",time_str(time,0),
            sat,iode);
        return NULL;
    }
    return nav->geph+j;
}
/* select sbas ephememeris ---------------------------------------------------*/
static seph_t *selseph(gtime_t time, int sat, const nav_t *nav)
{
    double t,tmax=MAXDTOE_SBS,tmin=tmax+1.0;
    int i,j=-1;

    trace(4,"selseph : time=%s sat=%2d\n",time_str(time,3),sat);

    for (i=0;i<nav->ns;i++) {
        if (nav->seph[i].sat!=sat) continue;
        if ((t=fabs(timediff(nav->seph[i].t0,time)))>tmax) continue;
        if (t<=tmin) {j=i; tmin=t;} /* toe closest to time */
    }
    if (j<0) {
        trace(3,"no sbas ephemeris     : %s sat=%2d\n",time_str(time,0),sat);
        return NULL;
    }
    return nav->seph+j;
}
/* satellite clock with broadcast ephemeris ----------------------------------*/
extern int ephclk(gtime_t time, gtime_t teph, int sat, const nav_t *nav,
    double *dts)
{
    eph_t  *eph;
    geph_t *geph;
    seph_t *seph;
    int sys;

    trace(4,"ephclk  : time=%s sat=%2d\n",time_str(time,3),sat);

    sys=satsys(sat,NULL);

    if (sys==SYS_GPS||sys==SYS_GAL||sys==SYS_QZS||sys==SYS_CMP||sys==SYS_IRN) {
        if (!(eph=seleph(teph,sat,-1,nav))) return 0;
        *dts=eph2clk(time,eph);
    }
    else if (sys==SYS_GLO) {
        if (!(geph=selgeph(teph,sat,-1,nav))) return 0;
        *dts=geph2clk(time,geph);
    }
    else if (sys==SYS_SBS) {
        if (!(seph=selseph(teph,sat,nav))) return 0;
        *dts=seph2clk(time,seph);
    }
    else return 0;

    return 1;
}
/* satellite position and clock by broadcast ephemeris -----------------------*/
extern int ephpos(gtime_t time, gtime_t teph, int sat, const nav_t *nav,
    int iode, double *rs, double *dts, double *var, int *svh)
{
    eph_t  *eph;
    geph_t *geph;
    seph_t *seph;
    double rst[3],dtst[1],tt=1E-3;
    int i,sys;

    trace(4,"ephpos  : time=%s sat=%2d iode=%d\n",time_str(time,3),sat,iode);

    sys=satsys(sat,NULL);

    *svh=-1;

    if (sys==SYS_GPS||sys==SYS_GAL||sys==SYS_QZS||sys==SYS_CMP||sys==SYS_IRN) {
        if (!(eph=seleph(teph,sat,iode,nav))) return 0;
        eph2pos(time,eph,rs,dts,var);
        time=timeadd(time,tt);
        eph2pos(time,eph,rst,dtst,var);
        *svh=eph->svh;
    }
    else if (sys==SYS_GLO) {
        if (!(geph=selgeph(teph,sat,iode,nav))) return 0;
        geph2pos(time,geph,rs,dts,var);
        time=timeadd(time,tt);
        geph2pos(time,geph,rst,dtst,var);
        *svh=geph->svh;
    }
    else if (sys==SYS_SBS) {
        if (!(seph=selseph(teph,sat,nav))) return 0;
        seph2pos(time,seph,rs,dts,var);
        time=timeadd(time,tt);
        seph2pos(time,seph,rst,dtst,var);
        *svh=seph->svh;
    }
    else return 0;

    /* satellite velocity and clock drift by differential approx */
    for (i=0;i<3;i++) rs[i+3]=(rst[i]-rs[i])/tt;
    dts[1]=(dtst[0]-dts[0])/tt;

    return 1;
}
/* get selected satellite ephemeris -------------------------------------------
* Get the selected satellite ephemeris.
* args   : int    sys       I   satellite system (SYS_???)
* return : selected ephemeris
*            refer setseleph()
*-----------------------------------------------------------------------------*/
extern int getseleph(int sys)
{
    switch (sys) {
    case SYS_GPS: return eph_sel[0];
    case SYS_GLO: return eph_sel[1];
    case SYS_GAL: return eph_sel[2];
    case SYS_QZS: return eph_sel[3];
    case SYS_CMP: return eph_sel[4];
    case SYS_IRN: return eph_sel[5];
    case SYS_SBS: return eph_sel[6];
    }
    return 0;
}