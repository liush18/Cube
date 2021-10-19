
#include "sat.h"

#define SQR(x)      ((x)*(x))
#define STD_BRDCCLK 30.0          /* error of broadcast clock (m) */

/* satellite position and clock ------------------------------------------------
* compute satellite position, velocity and clock
* args  :gtime_t time     I   time (gpst)
*          gtime_t teph     I   time to select ephemeris (gpst)
*          int    sat       I   satellite number
*          nav_t  *nav      I   navigation data
*          int    ephopt    I   ephemeris option (EPHOPT_???)
*          double *rs       O   sat position and velocity (ecef)
*                               {x,y,z,vx,vy,vz} (m|m/s)
*          double *dts      O   sat clock {bias,drift} (s|s/s)
*          double *var      O   sat position and clock error variance (m^2)
*          int    *svh      O   sat health flag (-1:correction not available)
* return:status (1:ok,0:error)
* notes :satellite position is referenced to antenna phase center
*          satellite clock does not include code bias correction (tgd or bgd)
*-----------------------------------------------------------------------------*/
static int satpos(gtime_t time, gtime_t teph, int sat, int ephopt,
    const nav_t *nav, double *rs, double *dts, double *var, int *svh)
{
    trace(4,"satpos: time=%s sat=%2d ephopt=%d\n",time_str(time,3),sat,ephopt);

    *svh=0;
    switch (ephopt) {
        case EPHOPT_BRDC: return ephpos(time,teph,sat,nav,-1,rs,dts,var,svh);
        case EPHOPT_PREC:
            if (!peph2pos(time,sat,nav,1,rs,dts,var)) break; else return 1;
    }
    *svh=-1;
    return 0;
}
/* satellite positions and clocks ----------------------------------------------
* compute satellite positions, velocities and clocks
* args  :gtime_t teph     I   time to select ephemeris (gpst)
*          obsd_t *obs      I   observation data
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation data
*          int    ephopt    I   ephemeris option (EPHOPT_???)
*          double *rs       O   satellite positions and velocities (ecef)
*          double *dts      O   satellite clocks
*          double *var      O   sat position and clock error variances (m^2)
*          int    *svh      O   sat health flag (-1:correction not available)
* return:none
* notes :rs [(0:2)+i*6]= obs[i] sat position {x,y,z} (m)
*          rs [(3:5)+i*6]= obs[i] sat velocity {vx,vy,vz} (m/s)
*          dts[(0:1)+i*2]= obs[i] sat clock {bias,drift} (s|s/s)
*          var[i]       =obs[i] sat position and clock error variance (m^2)
*          svh[i]       =obs[i] sat health flag
*          if no navigation data, set 0 to rs[], dts[], var[] and svh[]
*          satellite position and clock are values at signal transmission time
*          satellite position is referenced to antenna phase center
*          satellite clock does not include code bias correction (tgd or bgd)
*          any pseudorange and broadcast ephemeris are always needed to get
*          signal transmission time
*-----------------------------------------------------------------------------*/
extern void satposs(gtime_t teph, const obsd_t *obs, int n, const nav_t *nav,
    int ephopt, double *rs, double *dts, double *var, int *svh)
{
    gtime_t time[2*MAXOBS]={{0}};
    double dt,pr;
    int i,j;

    trace(3,"satposs:teph=%s n=%d ephopt=%d\n",time_str(teph,3),n,ephopt);

    for (i=0;i<n&&i<2*MAXOBS;i++) {
        for (j=0;j<6;j++) rs [j+i*6]=0.0;
        for (j=0;j<2;j++) dts[j+i*2]=0.0;
        var[i]=0.0; svh[i]=0;

        /* search any pseudorange */
        for (j=0,pr=0.0;j<NFREQ;j++) if ((pr=obs[i].P[j])!=0.0) break;

        if (j>=NFREQ) {
            trace(3,"no pseudorange %s sat=%2d\n",time_str(obs[i].time,3),obs[i].sat);
            continue;
        }
        /* transmission time by satellite clock */
        time[i]=timeadd(obs[i].time,-pr/CLIGHT);

        /* satellite clock bias by broadcast ephemeris */
        if (!ephclk(time[i],teph,obs[i].sat,nav,&dt)) {
            trace(3,"no broadcast clock %s sat=%2d\n",time_str(time[i],3),obs[i].sat);
            continue;
        }
        time[i]=timeadd(time[i],-dt);

        /* satellite position and clock at transmission time */
        if (!satpos(time[i],teph,obs[i].sat,ephopt,nav,rs+i*6,dts+i*2,var+i,
            svh+i)) {
            trace(3,"no ephemeris %s sat=%2d\n",time_str(time[i],3),obs[i].sat);
            continue;
        }
        /* if no precise clock available, use broadcast clock instead */
        if (dts[i*2]==0.0) {
            if (!ephclk(time[i],teph,obs[i].sat,nav,dts+i*2)) continue;
            dts[1+i*2]=0.0;
            *var=SQR(STD_BRDCCLK);
        }
    }
    for (i=0;i<n&&i<2*MAXOBS;i++) {
        trace(4,"%s sat=%2d rs=%13.3f %13.3f %13.3f dts=%12.3f var=%7.3f svh=%02X\n",
            time_str(time[i],6),obs[i].sat,rs[i*6],rs[1+i*6],rs[2+i*6],
            dts[i*2]*1E9,var[i],svh[i]);
    }
}

/* satellite positions and clocks ----------------------------------------------
* compute satellite positions, velocities and clocks for decoupled model
* args  :gtime_t teph     I   time to select ephemeris (gpst)
*          obsd_t *obs      I   observation data
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation data
*          double *rs       O   satellite positions and velocities (ecef)
*          double *dts      O   satellite clocks
*          double *vars     O   sat position error variances (m^2)
*          double *varc     O   sat clock error variances (m^2)
*          int    *svh      O   sat health flag (-1:correction not available)
* return:none
* notes :rs [(0:2)+i*6]=obs[i] sat position {x,y,z} (m)
*          rs [(3:5)+i*6]=obs[i] sat velocity {vx,vy,vz} (m/s)
*          dts[(0:2)+i*3]=obs[i] sat clock {p_bias,l_bias,m_bias} (s|m)
*          vars[i]       =obs[i] sat position error variance (m^2)
*          varc[(0:2)+i*3]= obs[i] sat clock error variance (m^2)
*          svh[i]        =obs[i] sat health flag
*-----------------------------------------------------------------------------*/
extern void sat_poss(gtime_t teph, char *rcv, const obsd_t *obs, int n, 
    const nav_t *nav, double *rs, double *dts, double *vars, double *varc, 
    int *svh)
{
    gtime_t time[2*MAXOBS]={{0}};
    double dt,pr;
    char id[8],str[32];
    int i,j;

    time2str(obs[0].time,str,2);

    for (i=0;i<n&&i<2*MAXOBS;i++) {
        for (j=0;j<6;j++) rs[j+i*6]=0.0;
        if (dts)  for (j=0;j<3;j++) dts[j+i*3]=0.0;
        if (varc) for (j=0;j<3;j++) varc[j+i*3]=0.0;
        vars[i]=0.0; svh[i]=0;

        /* search any pseudorange */
        for(j=0,pr=0.0;j<NFREQ;j++) if((pr=obs[i].P[j])!=0.0) break;

        satno2id(obs[i].sat,id);
        if(j>=NFREQ){
            trace(3,"%s rcv(%s) sat(%s) no pseudorange\n",str,rcv,id);
            continue;
        }
        /* transmission time by satellite clock */
        time[i]=timeadd(obs[i].time,-pr/CLIGHT);

        /* satellite clock bias by broadcast ephemeris */
        if (!ephclk(time[i],teph,obs[i].sat,nav,&dt)){
            trace(3,"%s rcv(%s) sat(%s) no broadcast clock\n",str,rcv,id);
            continue;
        }
        time[i]=timeadd(time[i],-dt);

        if (dts&&varc) {
            if (!peph_to_pos(time[i],obs[i].sat,nav,1,rs+i*6,dts+i*3,vars+i,varc+i*3))
                svh[i]=-1;
        }
        else {
            if (!peph_to_pos(time[i],obs[i].sat,nav,1,rs+i*6,NULL,vars+i,NULL))
                svh[i]=-1;
        }
    }

    if (dts&&varc) for (i=0;i<n&&i<2*MAXOBS;i++) {
        trace(4,"%s rcv=%s sat=%2d rs=%13.3f %13.3f %13.3f vars=%7.3f \
            dts=%12.3fns %12.3fns %12.3fm varc=%7.3f %7.3f %7.3f svh=%02X\n",
            str,rcv,obs[i].sat,rs[i*6],rs[1+i*6],rs[2+i*6],vars[i],
            dts[i*3]*1E9,dts[1+i*3]*1E9,dts[2+i*3],
            varc[i*3],varc[1+i*3],varc[2+i*3],svh[i]);
    }
    else for (i=0;i<n&&i<2*MAXOBS;i++) {
        trace(4,"%s rcv=%s sat=%2d rs=%13.3f %13.3f %13.3f vars=%7.3f svh=%02X\n",
            str,rcv,obs[i].sat,rs[i*6],rs[1+i*6],rs[2+i*6],vars[i],svh[i]);
    }
} 