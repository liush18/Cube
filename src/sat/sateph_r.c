/* GLO, ephemeris to satellite clock and position */
#include "sat.h"

#define SQR(x)   ((x)*(x))

#define RE_GLO      6378136.0       /* radius of earth (m)            ref [2] */
#define MU_GLO      3.9860044E14    /* gravitational constant         ref [2] */
#define J2_GLO      1.0826257E-3    /* 2nd zonal harmonic of geopot   ref [2] */
#define OMGE_GLO    7.292115E-5     /* earth angular velocity (rad/s) ref [2] */
#define ERREPH_GLO  5.0             /* error of glonass ephemeris (m) */
#define TSTEP       60.0            /* integration step glonass ephemeris (s) */

/* glonass orbit differential equations --------------------------------------*/
static void deq(const double *x, double *xdot, const double *acc)
{
    double a,b,c,r2=dot(x,x,3),r3=r2*sqrt(r2),omg2=SQR(OMGE_GLO);

    if (r2<=0.0) {
        xdot[0]=xdot[1]=xdot[2]=xdot[3]=xdot[4]=xdot[5]=0.0;
        return;
    }
    /* ref [2] A.3.1.2 with bug fix for xdot[4],xdot[5] */
    a=1.5*J2_GLO*MU_GLO*SQR(RE_GLO)/r2/r3; /* 3/2*J2*mu*Ae^2/r^5 */
    b=5.0*x[2]*x[2]/r2;                    /* 5*z^2/r^2 */
    c=-MU_GLO/r3-a*(1.0-b);                /* -mu/r^3-a(1-b) */
    xdot[0]=x[3]; xdot[1]=x[4]; xdot[2]=x[5];
    xdot[3]=(c+omg2)*x[0]+2.0*OMGE_GLO*x[4]+acc[0];
    xdot[4]=(c+omg2)*x[1]-2.0*OMGE_GLO*x[3]+acc[1];
    xdot[5]=(c-2.0*a)*x[2]+acc[2];
}
/* glonass position and velocity by numerical integration --------------------*/
static void glorbit(double t, double *x, const double *acc)
{
    double k1[6],k2[6],k3[6],k4[6],w[6];
    int i;

    deq(x,k1,acc); for (i=0;i<6;i++) w[i]=x[i]+k1[i]*t/2.0;
    deq(w,k2,acc); for (i=0;i<6;i++) w[i]=x[i]+k2[i]*t/2.0;
    deq(w,k3,acc); for (i=0;i<6;i++) w[i]=x[i]+k3[i]*t;
    deq(w,k4,acc);
    for (i=0;i<6;i++) x[i]+=(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i])*t/6.0;
}
/* glonass ephemeris to satellite clock bias -----------------------------------
* compute satellite clock bias with glonass ephemeris
* args   : gtime_t time     I   time by satellite clock (gpst)
*          geph_t *geph     I   glonass ephemeris
* return : satellite clock bias (s)
* notes  : see ref [2]
*-----------------------------------------------------------------------------*/
extern double geph2clk(gtime_t time, const geph_t *geph)
{
    double t,ts;
    int i;

    trace(4,"geph2clk: time=%s sat=%2d\n",time_str(time,3),geph->sat);

    t=ts=timediff(time,geph->toe);

    for (i=0;i<2;i++) {
        t=ts-(-geph->taun+geph->gamn*t);
    }
    return -geph->taun+geph->gamn*t;
}
/* glonass ephemeris to satellite position and clock bias ----------------------
* compute satellite position and clock bias with glonass ephemeris
* args  :gtime_t time     I   time (gpst)
*          geph_t *geph     I   glonass ephemeris
*          double *rs       O   satellite position {x,y,z} (ecef) (m)
*          double *dts      O   satellite clock bias (s)
*          double *var      O   satellite position and clock variance (m^2)
* return:none
* notes :see ref [2]
*-----------------------------------------------------------------------------*/
extern void geph2pos(gtime_t time, const geph_t *geph, double *rs, double *dts,
    double *var)
{
    double t,tt,x[6];
    int i;

    trace(4,"geph2pos: time=%s sat=%2d\n",time_str(time,3),geph->sat);

    t=timediff(time,geph->toe);

    *dts=-geph->taun+geph->gamn*t;

    for (i=0;i<3;i++) {
        x[i  ]=geph->pos[i];
        x[i+3]=geph->vel[i];
    }
    for (tt=t<0.0?-TSTEP:TSTEP;fabs(t)>1E-9;t-=tt) {
        if (fabs(t)<TSTEP) tt=t;
        glorbit(tt,x,geph->acc);
    }
    for (i=0;i<3;i++) rs[i]=x[i];

    *var=SQR(ERREPH_GLO);
}