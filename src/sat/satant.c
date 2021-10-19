
#include "sat.h"

#define SQR(x)      ((x)*(x))

/* satellite antenna phase center offset ---------------------------------------
* compute satellite antenna phase center offset in ecef
* args   : gtime_t time       I   time (gpst)
*          double *rs         I   satellite position and velocity (ecef)
*                                 {x,y,z,vx,vy,vz} (m|m/s)
*          int    sat         I   satellite number
*          nav_t  *nav        I   navigation data
*          double *dant       I   satellite antenna phase center offset (ecef)
*                                 {dx,dy,dz} (m) (iono-free LC value)
* return : none
* notes  : iono-free LC frequencies defined as follows:
*            GPS/QZSS : L1-L2
*            GLONASS  : G1-G2
*            Galileo  : E1-E5b
*            BDS      : B1I-B2I
*            NavIC    : L5-S
*-----------------------------------------------------------------------------*/
extern void satantoff(gtime_t time, const double *rs, int sat, const nav_t *nav,
    double *dant)
{
    const pcv_t *pcv=nav->pcvs+sat-1;
    double ex[3],ey[3],ez[3],es[3],r[3],rsun[3],gmst,erpv[5]={0},freq[2];
    double C1,C2,dant1,dant2;
    int i,sys;

    trace(4,"satantoff: time=%s sat=%2d\n",time_str(time,3),sat);

    dant[0]=dant[1]=dant[2]=0.0;

    /* sun position in ecef */
    sunmoonpos(gpst2utc(time),erpv,rsun,NULL,&gmst);

    /* unit vectors of satellite fixed coordinates */
    for (i=0;i<3;i++) r[i]=-rs[i];
    if (!normv3(r,ez)) return;
    for (i=0;i<3;i++) r[i]=rsun[i]-rs[i];
    if (!normv3(r,es)) return;
    cross3(ez,es,r);
    if (!normv3(r,ey)) return;
    cross3(ey,ez,ex);

    /* iono-free LC coefficients */
    sys=satsys(sat,NULL);
    if (sys==SYS_GPS||sys==SYS_QZS) { /* L1-L2 */
        freq[0]=FREQ1;
        freq[1]=FREQ2;
    }
    else if (sys==SYS_GLO) { /* G1-G2 */
        freq[0]=sat2freq(sat,CODE_L1C,nav);
        freq[1]=sat2freq(sat,CODE_L2C,nav);
    }
    else if (sys==SYS_GAL) { /* E1-E5b */
        freq[0]=FREQ1;
        freq[1]=FREQ7;
    }
    else if (sys==SYS_CMP) { /* B1I-B2I */
        freq[0]=FREQ1_CMP;
        freq[1]=FREQ2_CMP;
    }
    else if (sys==SYS_IRN) { /* B1I-B2I */
        freq[0]=FREQ5;
        freq[1]=FREQ9;
    }
    else return;

    C1= SQR(freq[0])/(SQR(freq[0])-SQR(freq[1]));
    C2=-SQR(freq[1])/(SQR(freq[0])-SQR(freq[1]));

    /* iono-free LC */
    for (i=0;i<3;i++) {
        dant1=pcv->off[0][0]*ex[i]+pcv->off[0][1]*ey[i]+pcv->off[0][2]*ez[i];
        dant2=pcv->off[1][0]*ex[i]+pcv->off[1][1]*ey[i]+pcv->off[1][2]*ez[i];
        dant[i]=C1*dant1+C2*dant2;
    }
}