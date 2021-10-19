
#include "sat.h"

#define SQR(x)      ((x)*(x))

#define NMAX        10              /* order of polynomial interpolation */
#define MAXDTE      900.0           /* max time difference to ephem time (s) */
#define EXTERR_CLK  1E-3            /* extrapolation error for clock (m/s) */
#define EXTERR_EPH  5E-7            /* extrapolation error for ephem (m/s^2) */

/* polynomial interpolation by Neville's algorithm ---------------------------*/
static double interppol(const double *x, double *y, int n)
{
    int i,j;

    for (j=1;j<n;j++) {
        for (i=0;i<n-j;i++) {
            y[i]=(x[i+j]*y[i]-x[i]*y[i+1])/(x[i+j]-x[i]);
        }
    }
    return y[0];
}
/* satellite position by precise ephemeris -----------------------------------*/
static int pephpos(gtime_t time, int sat, const nav_t *nav, double *rs,
    double *dts, double *vare, double *varc)
{
    double t[NMAX+1],p[3][NMAX+1],c[2],*pos,std=0.0,s[3],sinl,cosl;
    int i,j,k,index;

    trace(4,"pephpos:time=%s sat=%2d\n",time_str(time,3),sat);

    rs[0]=rs[1]=rs[2]=dts[0]=0.0;

    if (nav->ne<NMAX+1||
        timediff(time,nav->peph[0].time)<-MAXDTE||
        timediff(time,nav->peph[nav->ne-1].time)>MAXDTE) {
        trace(3,"no prec ephem %s sat=%2d\n",time_str(time,0),sat);
        return 0;
    }
    /* binary search */
    for (i=0,j=nav->ne-1;i<j;) {
        k=(i+j)/2;
        if (timediff(nav->peph[k].time,time)<0.0) i=k+1; else j=k;
    }
    index=i<=0?0:i-1;

    /* polynomial interpolation for orbit */
    i=index-(NMAX+1)/2;
    if (i<0) i=0; else if (i+NMAX>=nav->ne) i=nav->ne-NMAX-1;

    for (j=0;j<=NMAX;j++) {
        t[j]=timediff(nav->peph[i+j].time,time);
        if (norm(nav->peph[i+j].pos[sat-1],3)<=0.0) {
            trace(3,"prec ephem outage %s sat=%2d\n",time_str(time,0),sat);
            return 0;
        }
    }
    for (j=0;j<=NMAX;j++) {
        pos=nav->peph[i+j].pos[sat-1];
        /* correciton for earh rotation ver.2.4.0 */
        sinl=sin(OMGE*t[j]);
        cosl=cos(OMGE*t[j]);
        p[0][j]=cosl*pos[0]-sinl*pos[1];
        p[1][j]=sinl*pos[0]+cosl*pos[1];
        p[2][j]=pos[2];
    }
    for (i=0;i<3;i++) {
        rs[i]=interppol(t,p[i],NMAX+1);
    }
    if (vare) {
        for (i=0;i<3;i++) s[i]=nav->peph[index].std[sat-1][i];
        std=norm(s,3);

        /* extrapolation error for orbit */
        if      (t[0   ]>0.0) std+=EXTERR_EPH*SQR(t[0   ])/2.0;
        else if (t[NMAX]<0.0) std+=EXTERR_EPH*SQR(t[NMAX])/2.0;
        *vare=SQR(std);
    }
    /* linear interpolation for clock */
    t[0]=timediff(time,nav->peph[index  ].time);
    t[1]=timediff(time,nav->peph[index+1].time);
    c[0]=nav->peph[index  ].pos[sat-1][3];
    c[1]=nav->peph[index+1].pos[sat-1][3];

    if (t[0]<=0.0) {
        if ((dts[0]=c[0])!=0.0) {
            std=nav->peph[index].std[sat-1][3]*CLIGHT-EXTERR_CLK*t[0];
        }
    }
    else if (t[1]>=0.0) {
        if ((dts[0]=c[1])!=0.0) {
            std=nav->peph[index+1].std[sat-1][3]*CLIGHT+EXTERR_CLK*t[1];
        }
    }
    else if (c[0]!=0.0&&c[1]!=0.0) {
        dts[0]=(c[1]*t[0]-c[0]*t[1])/(t[0]-t[1]);
        i=t[0]<-t[1]?0:1;
        std=nav->peph[index+i].std[sat-1][3]+EXTERR_CLK*fabs(t[i]);
    }
    else {
        dts[0]=0.0;
    }
    if (varc) *varc=SQR(std);
    return 1;
}
/* satellite clock by precise clock ------------------------------------------*/
static int pephclk(gtime_t time, int sat, const nav_t *nav, double *dts,
    double *varc)
{
    double t[2],c[2],std;
    int i,j,k,index;

    trace(4,"pephclk:time=%s sat=%2d\n",time_str(time,3),sat);

    if (nav->nc<2||
        timediff(time,nav->pclk[0].time)<-MAXDTE||
        timediff(time,nav->pclk[nav->nc-1].time)>MAXDTE) {
        trace(3,"no prec clock %s sat=%2d\n",time_str(time,0),sat);
        return 1;
    }
    /* binary search */
    for (i=0,j=nav->nc-1;i<j;) {
        k=(i+j)/2;
        if (timediff(nav->pclk[k].time,time)<0.0) i=k+1; else j=k;
    }
    index=i<=0?0:i-1;

    /* linear interpolation for clock */
    t[0]=timediff(time,nav->pclk[index  ].time);
    t[1]=timediff(time,nav->pclk[index+1].time);
    c[0]=nav->pclk[index  ].clk[sat-1][0];
    c[1]=nav->pclk[index+1].clk[sat-1][0];

    /* t[0]<=0.0, EXTERR_CLK*t[0]>=0.0 */
    if (t[0]<=0.0) {
        if ((dts[0]=c[0])==0.0) return 0;
        std=nav->pclk[index].std[sat-1][0]*CLIGHT-EXTERR_CLK*t[0];
    }
    else if (t[1]>=0.0) {
        if ((dts[0]=c[1])==0.0) return 0;
        std=nav->pclk[index+1].std[sat-1][0]*CLIGHT+EXTERR_CLK*t[1];
    }
    else if (c[0]!=0.0&&c[1]!=0.0) {
        dts[0]=(c[1]*t[0]-c[0]*t[1])/(t[0]-t[1]);
        i=t[0]<-t[1]?0:1;
        std=nav->pclk[index+i].std[sat-1][0]*CLIGHT+EXTERR_CLK*fabs(t[i]);
    }
    else {
        trace(3,"prec clock outage %s sat=%2d\n",time_str(time,0),sat);
        return 0;
    }
    if (varc) *varc=SQR(std);
    return 1;
}
/* satellite position/clock by precise ephemeris/clock -------------------------
* compute satellite position/clock with precise ephemeris/clock
* args  :gtime_t time       I   time (gpst)
*          int    sat         I   satellite number
*          nav_t  *nav        I   navigation data
*          int    opt         I   sat postion option
*                                 (0: center of mass, 1: antenna phase center)
*          double *rs         O   sat position and velocity (ecef)
*                                 {x,y,z,vx,vy,vz} (m|m/s)
*          double *dts        O   sat clock {bias,drift} (s|s/s)
*          double *var        IO  sat position and clock error variance (m)
*                                 (NULL: no output)
* return:status (1:ok,0:error or data outage)
* notes :clock includes relativistic correction but does not contain code bias
*          before calling the function, nav->peph, nav->ne, nav->pclk and
*          nav->nc must be set by calling readsp3(), readrnx() or readrnxt()
*          if precise clocks are not set, clocks in sp3 are used instead
*-----------------------------------------------------------------------------*/
extern int peph2pos(gtime_t time, int sat, const nav_t *nav, int opt,
    double *rs, double *dts, double *var)
{
    gtime_t time_tt;
    double rss[3],rst[3],dtss[1],dtst[1],dant[3]={0},vare=0.0,varc=0.0,tt=1E-3;
    int i;

    trace(4,"peph2pos: time=%s sat=%2d opt=%d\n",time_str(time,3),sat,opt);

    if (sat<=0||MAXSAT<sat) return 0;

    /* satellite position and clock bias */
    if (!pephpos(time,sat,nav,rss,dtss,&vare,&varc)||
        !pephclk(time,sat,nav,dtss,&varc)) return 0;

    time_tt=timeadd(time,tt);
    if (!pephpos(time_tt,sat,nav,rst,dtst,NULL,NULL)||
        !pephclk(time_tt,sat,nav,dtst,NULL)) return 0;

    /* satellite antenna offset correction */
    if (opt) {
        satantoff(time,rss,sat,nav,dant);
    }
    for (i=0;i<3;i++) {
        rs[i  ]=rss[i]+dant[i];
        rs[i+3]=(rst[i]-rss[i])/tt;
    }
    /* relativistic effect correction (remove to relativity()) */
    if (dtss[0]!=0.0) {
        dts[0]=dtss[0];//-2.0*dot(rs,rs+3,3)/CLIGHT/CLIGHT;
        dts[1]=(dtst[0]-dtss[0])/tt;
    }
    else { /* no precise clock */
        dts[0]=dts[1]=0.0;
    }
    if (var) *var=vare+varc;

    return 1;
}

/* precise empheris to position and clock, for decoupled clock model
*only position needed for decoupled clock estimated
*position and three clocks for every satellite needed for decoupled clock ppp
 *----------------------------------------------------------------------------*/
/* satellite position by precise ephemeris -----------------------------------*/
static int peph_pos(gtime_t time, int sat, const nav_t *nav, double *rs,
    double *vare)
{
    double t[NMAX+1],p[3][NMAX+1],*pos,std=0.0,s[3],sinl,cosl;
    int i,j,k,index;

    trace(4,"peph_pos: time=%s sat=%2d\n",time_str(time,3),sat);

    rs[0]=rs[1]=rs[2]=0.0;

    if (nav->ne<NMAX+1||
        timediff(time,nav->peph[0].time)<-MAXDTE||
        timediff(time,nav->peph[nav->ne-1].time)>MAXDTE) {
        trace(3,"no prec ephem %s sat=%2d\n",time_str(time,0),sat);
        return 0;
    }
    /* binary search */
    for (i=0,j=nav->ne-1;i<j;) {
        k=(i+j)/2;
        if (timediff(nav->peph[k].time,time)<0.0) i=k+1; else j=k;
    }
    index=i<=0?0:i-1;

    /* polynomial interpolation for orbit */
    i=index-(NMAX+1)/2;
    if(i<0) i=0; else if (i+NMAX>=nav->ne) i=nav->ne-NMAX-1;

    for (j=0;j<=NMAX;j++) {
        t[j]=timediff(nav->peph[i+j].time,time);
        if (norm(nav->peph[i+j].pos[sat-1],3)<=0.0) {
            trace(3,"prec ephem outage %s sat=%2d\n",time_str(time,0),sat);
            return 0;
        }
    }
    for (j=0;j<=NMAX;j++) {
        pos=nav->peph[i+j].pos[sat-1];

        /* correciton for earh rotation ver.2.4.0 */
        sinl=sin(OMGE*t[j]);
        cosl=cos(OMGE*t[j]);
        p[0][j]=cosl*pos[0]-sinl*pos[1];
        p[1][j]=sinl*pos[0]+cosl*pos[1];
        p[2][j]=pos[2];
    }
    for (i=0;i<3;i++) {
        rs[i]=interppol(t,p[i],NMAX+1);
    }
    if (vare) {
        for (i=0;i<3;i++) s[i]=nav->peph[index].std[sat-1][i];
        std=norm(s,3);

        /* extrapolation error for orbit */
        if      (t[0   ]>0.0) std+=EXTERR_EPH*SQR(t[0   ])/2.0;
        else if (t[NMAX]<0.0) std+=EXTERR_EPH*SQR(t[NMAX])/2.0;
        *vare=SQR(std);
    }
    return 1;
}
/* satellite clock by precise clock ------------------------------------------*/
static int peph_clk(gtime_t time, int sat, const nav_t *nav, double *dts,
    double *varc)
{
    double t[2],c[2],std,fact;
    int i,j,k,index;

    trace(4,"pephclk:time=%s sat=%2d\n",time_str(time,3),sat);

    if (nav->nc<2||
        timediff(time,nav->pclk[0].time)<-MAXDTE||
        timediff(time,nav->pclk[nav->nc-1].time)>MAXDTE) {
        trace(2,"no prec clock %s sat=%2d\n",time_str(time,0),sat);
        return 0;
    }
    /* binary search */
    for (i=0,j=nav->nc-1;i<j;) {
        k=(i+j)/2;
        if (timediff(nav->pclk[k].time,time)<0.0) i=k+1; else j=k;
    }
    index=i<=0?0:i-1;

    /* linear interpolation for clock */
    t[0]=timediff(time,nav->pclk[index  ].time);
    t[1]=timediff(time,nav->pclk[index+1].time);

    /* three satellite clocks */
    for (i=0;i<3;i++) {
        c[0]=nav->pclk[index  ].clk[sat-1][i];
        c[1]=nav->pclk[index+1].clk[sat-1][i];

        if (i<2) fact=CLIGHT;
        else fact=LAM4;
        if (t[0]<=0.0) {
            if ((dts[i]=c[0])==0.0) return 0;
            std=nav->pclk[index].std[sat-1][i]*fact-EXTERR_CLK*t[0];
        }
        else if (t[1]>=0.0) {
            if ((dts[i]=c[1])==0.0) return 0;
            std=nav->pclk[index+1].std[sat-1][i]*fact+EXTERR_CLK*t[1];
        }
        else if (c[0]!=0.0&&c[1]!=0.0) {
            dts[i]=(c[1]*t[0]-c[0]*t[1])/(t[0]-t[1]);
            k=t[0]<-t[1]?0:1;
            std=nav->pclk[index+k].std[sat-1][i]*fact+EXTERR_CLK*fabs(t[k]);
        }
        else {
            trace(2,"prec clock outage %s sat=%2d\n",time_str(time,0),sat);
            return 0;
        }
        if (varc) varc[i]=SQR(std);
    }
    return 1;
}
/* satellite position by precise ephemeris -----------------------------------*/
extern int peph_to_pos(gtime_t time, int sat, const nav_t *nav, int opt,
    double *rs, double *dts, double *vars, double *varc)
{
    gtime_t time_tt;
    double rss[3],rst[3],dtss[3],dant[3]={0},tt=1E-3;
    int i;

    trace(4,"peph2pos: time=%s sat=%2d opt=%d\n",time_str(time,3),sat,opt);

    if (sat<=0||MAXSAT<sat) return 0;

    time_tt=timeadd(time,tt);

    /* satellite position and clock bias */
    if (!peph_pos(time,sat,nav,rss,vars)) return 0;
    if (!peph_pos(time_tt,sat,nav,rst,NULL)) return 0;

    /* satellite antenna offset correction */
    if (opt) satantoff(time,rss,sat,nav,dant);

    for (i=0;i<3;i++) {
        rs[i  ]=rss[i]+dant[i];
        rs[i+3]=(rst[i]-rss[i])/tt;
    }

    if (dts) {
        if (!peph_clk(time,sat,nav,dtss,varc)) return 0;

        if (dtss[0]!=0.0&&dtss[1]!=0.0&&dtss[2]!=0.0) {
            /* corrected in relativity() */
            for (i=0;i<3;i++) {
                dts[i]=dtss[i];
                //dts[i]-=i<2?2.0*dot(rs,rs+3,3)/CLIGHT/CLIGHT:0.0;
            }
        }
        else {
            dts[0]=dts[1]=dts[2]=0.0;
            return 0;
        } 
    }

    return 1;
}