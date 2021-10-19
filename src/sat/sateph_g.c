/* -----------------------------------------------------------------------------
* GPS, GAL and BDS, ephemeris to satellite clock and position
* references :
*     [1] IS-GPS-200K, Navstar GPS Space Segment/Navigation User Interfaces,
*         May 6, 2019
*     [2] Global Navigation Satellite System GLONASS, Interface Control Document
*         Navigational radiosignal In bands L1, L2, (Version 5.1), 2008
*     [3] RTCA/DO-229C, Minimum operational performance standards for global
*         positioning system/wide area augmentation system airborne equipment,
*         RTCA inc, November 28, 2001
*     [4] RTCM Paper, April 12, 2010, Proposed SSR Messages for SV Orbit Clock,
*         Code Biases, URA
*     [5] RTCM Paper 012-2009-SC104-528, January 28, 2009 (previous ver of [4])
*     [6] RTCM Paper 012-2009-SC104-582, February 2, 2010 (previous ver of [4])
*     [7] European GNSS (Galileo) Open Service Signal In Space Interface Control
*         Document, Issue 1.3, December, 2016
*     [8] Quasi-Zenith Satellite System Interface Specification Satellite
*         Positioning, Navigation and Timing Service (IS-QZSS-PNT-003), Cabinet
*         Office, November 5, 2018
*     [9] BeiDou navigation satellite system signal in space interface control
*         document open service signal B1I (version 3.0), China Satellite
*         Navigation office, February, 2019
*     [10] RTCM Standard 10403.3, Differential GNSS (Global Navigation
*         Satellite Systems) Services - version 3, October 7, 2016
*-----------------------------------------------------------------------------*/
#include "sat.h"

#define SQR(x)   ((x)*(x))

#define MU_GPS   3.9860050E14     /* gravitational constant         ref [1] */
#define MU_GAL   3.986004418E14   /* earth gravitational constant   ref [7] */
#define MU_CMP   3.986004418E14   /* earth gravitational constant   ref [9] */

#define OMGE_GAL 7.2921151467E-5  /* earth angular velocity (rad/s) ref [7] */
#define OMGE_CMP 7.292115E-5      /* earth angular velocity (rad/s) ref [9] */

#define SIN_5 -0.0871557427476582 /* sin(-5.0 deg) */
#define COS_5  0.9961946980917456 /* cos(-5.0 deg) */

#define RTOL_KEPLER 1E-13         /* relative tolerance for Kepler equation */

#define DEFURASSR 0.15            /* default accurary of ssr corr (m) */
#define STD_GAL_NAPA 500.0        /* error of galileo ephemeris for NAPA (m) */

#define MAX_ITER_KEPLER 30        /* max number of iteration of Kelpler */

/* variance by ura ephemeris -------------------------------------------------*/
extern double var_uraeph(int sys, int ura)
{
    const double ura_value[]={   
        2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,
        3072.0,6144.0
    };
    if (sys==SYS_GAL) { /* galileo sisa (ref [7] 5.1.11) */
        if (ura<= 49) return SQR(ura*0.01);
        if (ura<= 74) return SQR(0.5+(ura- 50)*0.02);
        if (ura<= 99) return SQR(1.0+(ura- 75)*0.04);
        if (ura<=125) return SQR(2.0+(ura-100)*0.16);
        return SQR(STD_GAL_NAPA);
    }
    else { /* gps ura (ref [1] 20.3.3.3.1.1) */
        return ura<0||14<ura?SQR(6144.0):SQR(ura_value[ura]);
    }
}
/* variance by ura ssr (ref [10] table 3.3-1 DF389) --------------------------*/
static double var_urassr(int ura)
{
    double std;
    if (ura<= 0) return SQR(DEFURASSR);
    if (ura>=63) return SQR(5.4665);
    std=(pow(3.0,(ura>>3)&7)*(1.0+(ura&7)/4.0)-1.0)*1E-3;
    return SQR(std);
}
/* broadcast ephemeris to satellite clock bias ---------------------------------
* compute satellite clock bias with broadcast ephemeris (gps, galileo, qzss)
* args   : gtime_t time     I   time by satellite clock (gpst)
*          eph_t *eph       I   broadcast ephemeris
* return : satellite clock bias (s) without relativeity correction
* notes  : see ref [1],[7],[8]
*          satellite clock does not include relativity correction and tdg
*-----------------------------------------------------------------------------*/
extern double eph2clk(gtime_t time, const eph_t *eph)
{
    double t,ts;
    int i;

    trace(4,"eph2clk : time=%s sat=%2d\n",time_str(time,3),eph->sat);

    t=ts=timediff(time,eph->toc);

    for (i=0;i<2;i++) {
        t=ts-(eph->f0+eph->f1*t+eph->f2*t*t);
    }
    return eph->f0+eph->f1*t+eph->f2*t*t;
}
/* broadcast ephemeris to satellite position and clock bias --------------------
* compute satellite position and clock bias with broadcast ephemeris (gps,
* galileo, qzss)
* args   : gtime_t time     I   time (gpst)
*          eph_t *eph       I   broadcast ephemeris
*          double *rs       O   satellite position (ecef) {x,y,z} (m)
*          double *dts      O   satellite clock bias (s)
*          double *var      O   satellite position and clock variance (m^2)
* return : none
* notes  : see ref [1],[7],[8]
*          satellite clock includes relativity correction without code bias
*          (tgd or bgd)
*-----------------------------------------------------------------------------*/
extern void eph2pos(gtime_t time, const eph_t *eph, double *rs, double *dts,
    double *var)
{
    double tk,M,E,Ek,sinE,cosE,u,r,i,O,sin2u,cos2u,x,y,sinO,cosO,cosi,mu,omge;
    double xg,yg,zg,sino,coso;
    int n,sys,prn;

    trace(4,"eph2pos : time=%s sat=%2d\n",time_str(time,3),eph->sat);

    if (eph->A<=0.0) {
        rs[0]=rs[1]=rs[2]=*dts=*var=0.0;
        return;
    }
    tk=timediff(time,eph->toe);

    switch ((sys=satsys(eph->sat,&prn))) {
    case SYS_GAL: mu=MU_GAL; omge=OMGE_GAL; break;
    case SYS_CMP: mu=MU_CMP; omge=OMGE_CMP; break;
    default:      mu=MU_GPS; omge=OMGE;     break;
    }
    M=eph->M0+(sqrt(mu/(eph->A*eph->A*eph->A))+eph->deln)*tk;

    for (n=0,E=M,Ek=0.0;fabs(E-Ek)>RTOL_KEPLER&&n<MAX_ITER_KEPLER;n++) {
        Ek=E; E-=(E-eph->e*sin(E)-M)/(1.0-eph->e*cos(E));
    }
    if (n>=MAX_ITER_KEPLER) {
        trace(2,"eph2pos: kepler iteration overflow sat=%2d\n",eph->sat);
        return;
    }
    sinE=sin(E); cosE=cos(E);

    trace(4,"kepler: sat=%2d e=%8.5f n=%2d del=%10.3e\n",eph->sat,eph->e,n,E-Ek);

    u=atan2(sqrt(1.0-eph->e*eph->e)*sinE,cosE-eph->e)+eph->omg;
    r=eph->A*(1.0-eph->e*cosE);
    i=eph->i0+eph->idot*tk;
    sin2u=sin(2.0*u); cos2u=cos(2.0*u);
    u+=eph->cus*sin2u+eph->cuc*cos2u;
    r+=eph->crs*sin2u+eph->crc*cos2u;
    i+=eph->cis*sin2u+eph->cic*cos2u;
    x=r*cos(u); y=r*sin(u); cosi=cos(i);

    /* beidou geo satellite */
    if (sys==SYS_CMP&&(prn<=5||prn>=59)) { /* ref [9] table 4-1 */
        O=eph->OMG0+eph->OMGd*tk-omge*eph->toes;
        sinO=sin(O); cosO=cos(O);
        xg=x*cosO-y*cosi*sinO;
        yg=x*sinO+y*cosi*cosO;
        zg=y*sin(i);
        sino=sin(omge*tk); coso=cos(omge*tk);
        rs[0]= xg*coso+yg*sino*COS_5+zg*sino*SIN_5;
        rs[1]=-xg*sino+yg*coso*COS_5+zg*coso*SIN_5;
        rs[2]=-yg*SIN_5+zg*COS_5;
    }
    else {
        O=eph->OMG0+(eph->OMGd-omge)*tk-omge*eph->toes;
        sinO=sin(O); cosO=cos(O);
        rs[0]=x*cosO-y*cosi*sinO;
        rs[1]=x*sinO+y*cosi*cosO;
        rs[2]=y*sin(i);
    }
    tk=timediff(time,eph->toc);
    *dts=eph->f0+eph->f1*tk+eph->f2*tk*tk;

    /* relativity correction */
    *dts-=2.0*sqrt(mu*eph->A)*eph->e*sinE/SQR(CLIGHT);

    /* position and clock error variance */
    *var=var_uraeph(sys,eph->sva);
}