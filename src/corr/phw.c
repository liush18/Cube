
#include "corr.h"

/* nominal yaw-angle ---------------------------------------------------------*/
static double yaw_nominal(double beta, double mu)
{
    if (fabs(beta)<1E-12&&fabs(mu)<1E-12) return PI;
    return atan2(-tan(beta), sin(mu))+PI;
}
/* yaw-angle of satellite ----------------------------------------------------*/
extern int yaw_angle(int sat, const char* type, int opt, double beta, double mu,
    double* yaw)
{
    *yaw=yaw_nominal(beta, mu);
    return 1;
}
/* satellite attitude model --------------------------------------------------*/
static int sat_yaw(gtime_t time, int sat, const char* type, int opt,
    const double* rs, double* exs, double* eys)
{
    double rsun[3], ri[6], es[3], esun[3], n[3], p[3], en[3], ep[3], ex[3], E, beta, mu;
    double yaw, cosy, siny, erpv[5]={ 0 };
    int i;

    sunmoonpos(gpst2utc(time), erpv, rsun, NULL, NULL);

    /* beta and orbit angle */
    matcpy(ri, rs, 6, 1);
    ri[3]-=OMGE*ri[1];
    ri[4]+=OMGE*ri[0];
    cross3(ri, ri+3, n);
    cross3(rsun, n, p);
    if (!normv3(rs, es)||!normv3(rsun, esun)||!normv3(n, en) ||
        !normv3(p, ep)) return 0;
    beta=PI/2.0-acos(dot(esun, en, 3));
    E=acos(dot(es, ep, 3));
    mu=PI/2.0+(dot(es, esun, 3)<=0?-E:E);
    if (mu<-PI/2.0) mu+=2.0*PI;
    else if (mu>=PI/2.0) mu-=2.0*PI;

    /* yaw-angle of satellite */
    if (!yaw_angle(sat, type, opt, beta, mu, &yaw)) return 0;

    /* satellite fixed x,y-vector */
    cross3(en, es, ex);
    cosy=cos(yaw);
    siny=sin(yaw);
    for (i=0; i<3; i++) {
        exs[i]=-siny*en[i]+cosy*ex[i];
        eys[i]=-cosy*en[i]-siny*ex[i];
    }
    return 1;
}
/* phase windup model --------------------------------------------------------*/
extern int model_phw(gtime_t time, int sat, const char* type, int opt,
    const double* rs, const double* rr, double* phw)
{
    double exs[3], eys[3], ek[3], exr[3], eyr[3], eks[3], ekr[3], E[9];
    double dr[3], ds[3], drs[3], r[3], pos[3], cosp, ph;
    int i;

    if (opt<=0) return 1; /* no phase windup */

    /* satellite yaw attitude model */
    if (!sat_yaw(time, sat, type, opt, rs, exs, eys)) return 0;

    /* unit vector satellite to receiver */
    for (i=0; i<3; i++) r[i]=rr[i]-rs[i];
    if (!normv3(r, ek)) return 0;

    /* unit vectors of receiver antenna */
    ecef2pos(rr, pos);
    xyz2enu(pos, E);
    exr[0]=E[1]; exr[1]=E[4]; exr[2]=E[7]; /* x=north */
    eyr[0]=-E[0]; eyr[1]=-E[3]; eyr[2]=-E[6]; /* y=west  */

    /* phase windup effect */
    cross3(ek, eys, eks);
    cross3(ek, eyr, ekr);
    for (i=0; i<3; i++) {
        ds[i]=exs[i]-ek[i]*dot(ek, exs, 3)-eks[i];
        dr[i]=exr[i]-ek[i]*dot(ek, exr, 3)+ekr[i];
    }
    cosp=dot(ds, dr, 3)/norm(ds, 3)/norm(dr, 3);
    if (cosp<-1.0) cosp=-1.0;
    else if (cosp>1.0) cosp=1.0;
    ph=acos(cosp)/2.0/PI;
    cross3(ds, dr, drs);
    if (dot(ek, drs, 3)<0.0) ph=-ph;

    *phw=ph+floor(*phw-ph+0.5); /* in cycle */
    return 1;
}
/* phase windup correction -----------------------------------------------------
* phase windup correction (ref [7] 5.1.2)
* args  :gtime_t time     I   time (GPST)
*          double  *rs      I   satellite position (ecef) {x,y,z} (m)
*          double  *rr      I   receiver  position (ecef) {x,y,z} (m)
*          double  *phw     IO  phase windup correction (cycle)
* return:none
* notes :the previous value of phase windup correction should be set to *phw
*          as an input. the function assumes windup correction has no jump more
*          than 0.5 cycle.
*-----------------------------------------------------------------------------*/
extern void windupcorr(gtime_t time, const double *rs, const double *rr,
                       double *phw)
{
    double ek[3],exs[3],eys[3],ezs[3],ess[3],exr[3],eyr[3],eks[3],ekr[3],E[9];
    double dr[3],ds[3],drs[3],r[3],pos[3],rsun[3],cosp,ph,erpv[5]={0};
    int i;
    
    trace(4,"windupcorr: time=%s\n",time_str(time,0));
    
    /* sun position in ecef */
    sunmoonpos(gpst2utc(time),erpv,rsun,NULL,NULL);
    
    /* unit vector satellite to receiver */
    for (i=0;i<3;i++) r[i]=rr[i]-rs[i];
    if (!normv3(r,ek)) return;
    
    /* unit vectors of satellite antenna */
    for (i=0;i<3;i++) r[i]=-rs[i];
    if (!normv3(r,ezs)) return;
    for (i=0;i<3;i++) r[i]=rsun[i]-rs[i];
    if (!normv3(r,ess)) return;
    cross3(ezs,ess,r);
    if (!normv3(r,eys)) return;
    cross3(eys,ezs,exs);
    
    /* unit vectors of receiver antenna */
    ecef2pos(rr,pos);
    xyz2enu(pos,E);
    exr[0]= E[1]; exr[1]= E[4]; exr[2]= E[7]; /* x=north */
    eyr[0]=-E[0]; eyr[1]=-E[3]; eyr[2]=-E[6]; /* y=west  */
    
    /* phase windup effect */
    cross3(ek,eys,eks);
    cross3(ek,eyr,ekr);
    for (i=0;i<3;i++) {
        ds[i]=exs[i]-ek[i]*dot(ek,exs,3)-eks[i];
        dr[i]=exr[i]-ek[i]*dot(ek,exr,3)+ekr[i];
    }
    cosp=dot(ds,dr,3)/norm(ds,3)/norm(dr,3);
    if      (cosp<-1.0) cosp=-1.0;
    else if (cosp> 1.0) cosp= 1.0;
    ph=acos(cosp)/2.0/PI;
    cross3(ds,dr,drs);
    if (dot(ek,drs,3)<0.0) ph=-ph;
    
    *phw=ph+floor(*phw-ph+0.5); /* in cycle */
}