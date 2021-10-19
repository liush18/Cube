
#include "corr.h"

#define SQR(x)      ((x)*(x))
#define VAR_NOTEC   SQR(30.0)   /* variance of no tec */
#define MIN_EL      0.0         /* min elevation angle (rad) */
#define MIN_HGT     -1000.0     /* min user height (m) */

/* ionosphere model ------------------------------------------------------------
* compute ionospheric delay by broadcast ionosphere model (klobuchar model)
* args   : gtime_t t        I   time (gpst)
*          double *ion      I   iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3}
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
* return : ionospheric delay (L1) (m)
*-----------------------------------------------------------------------------*/
extern double ionmodel(gtime_t t, const double *ion, const double *pos,
                       const double *azel)
{
    const double ion_default[]={ /* 2004/1/1 */
        0.1118E-07,-0.7451E-08,-0.5961E-07, 0.1192E-06,
        0.1167E+06,-0.2294E+06,-0.1311E+06, 0.1049E+07
    };
    double tt,f,psi,phi,lam,amp,per,x;
    int week;
    
    if (pos[2]<-1E3||azel[1]<=0) return 0.0;
    if (norm(ion,8)<=0.0) ion=ion_default;
    
    /* earth centered angle (semi-circle) */
    psi=0.0137/(azel[1]/PI+0.11)-0.022;
    
    /* subionospheric latitude/longitude (semi-circle) */
    phi=pos[0]/PI+psi*cos(azel[0]);
    if      (phi> 0.416) phi= 0.416;
    else if (phi<-0.416) phi=-0.416;
    lam=pos[1]/PI+psi*sin(azel[0])/cos(phi*PI);
    
    /* geomagnetic latitude (semi-circle) */
    phi+=0.064*cos((lam-1.617)*PI);
    
    /* local time (s) */
    tt=43200.0*lam+time2gpst(t,&week);
    tt-=floor(tt/86400.0)*86400.0; /* 0<=tt<86400 */
    
    /* slant factor */
    f=1.0+16.0*pow(0.53-azel[1]/PI,3.0);
    
    /* ionospheric delay */
    amp=ion[0]+phi*(ion[1]+phi*(ion[2]+phi*ion[3]));
    per=ion[4]+phi*(ion[5]+phi*(ion[6]+phi*ion[7]));
    amp=amp<    0.0?    0.0:amp;
    per=per<72000.0?72000.0:per;
    x=2.0*PI*(tt-50400.0)/per;
    
    return CLIGHT*f*(fabs(x)<1.57?5E-9+amp*(1.0+x*x*(-0.5+x*x/24.0)):5E-9);
}
/* ionosphere mapping function -------------------------------------------------
* compute ionospheric delay mapping function by single layer model
* args   : double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
* return : ionospheric mapping function
*-----------------------------------------------------------------------------*/
extern double ionmapf(const double *pos, const double *azel)
{
    if (pos[2]>=HION) return 1.0;
    return 1.0/cos(asin((RE_WGS84+pos[2])/(RE_WGS84+HION)*sin(PI/2.0-azel[1])));
}
/* ionospheric pierce point position -------------------------------------------
* compute ionospheric pierce point (ipp) position and slant factor
* args   : double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          double re        I   earth radius (km)
*          double hion      I   altitude of ionosphere (km)
*          double *posp     O   pierce point position {lat,lon,h} (rad,m)
* return : slant factor
* notes  : see ref [2], only valid on the earth surface
*          fixing bug on ref [2] A.4.4.10.1 A-22,23
*-----------------------------------------------------------------------------*/
extern double ionppp(const double *pos, const double *azel, double re,
                     double hion, double *posp)
{
    double cosaz,rp,ap,sinap,tanap;
    
    rp=re/(re+hion)*cos(azel[1]);
    ap=PI/2.0-azel[1]-asin(rp);
    sinap=sin(ap);
    tanap=tan(ap);
    cosaz=cos(azel[0]);
    posp[0]=asin(sin(pos[0])*cos(ap)+cos(pos[0])*sinap*cosaz); /* latitude */
    /* longitude */
    if ((pos[0]> 70.0*D2R&& tanap*cosaz>tan(PI/2.0-pos[0]))||
        (pos[0]<-70.0*D2R&&-tanap*cosaz>tan(PI/2.0+pos[0]))) {
        posp[1]=pos[1]+PI-asin(sinap*sin(azel[0])/cos(posp[0]));
    }
    else {
        posp[1]=pos[1]+asin(sinap*sin(azel[0])/cos(posp[0]));
    }
    return 1.0/sqrt(1.0-rp*rp);
}
/* data index (i:lat,j:lon,k:hgt) --------------------------------------------*/
static int dataindex(int i, int j, int k, const int *ndata)
{
    if (i<0||ndata[0]<=i||j<0||ndata[1]<=j||k<0||ndata[2]<=k) return -1;
    return i+ndata[0]*(j+ndata[1]*k);
}
/* interpolate tec grid data -------------------------------------------------*/
static int interptec(const tec_t *tec, int k, const double *posp, double *value,
    double *rms)
{
    double dlat,dlon,a,b,d[4]={0},r[4]={0};
    int i,j,n,index;

    trace(3,"interptec: k=%d posp=%.2f %.2f\n",k,posp[0]*R2D,posp[1]*R2D);
    *value=*rms=0.0;

    if (tec->lats[2]==0.0||tec->lons[2]==0.0) return 0;

    dlat=posp[0]*R2D-tec->lats[0];
    dlon=posp[1]*R2D-tec->lons[0];
    if (tec->lons[2]>0.0) dlon-=floor( dlon/360)*360.0; /*  0<=dlon<360 */
    else                  dlon+=floor(-dlon/360)*360.0; /* -360<dlon<=0 */

    a=dlat/tec->lats[2];
    b=dlon/tec->lons[2];
    i=(int)floor(a); a-=i;
    j=(int)floor(b); b-=j;

    /* get gridded tec data */
    for (n=0;n<4;n++) {
        if ((index=dataindex(i+(n%2),j+(n<2?0:1),k,tec->ndata))<0) continue;
        d[n]=tec->data[index];
        r[n]=tec->rms [index];
    }
    if (d[0]>0.0&&d[1]>0.0&&d[2]>0.0&&d[3]>0.0) {

        /* bilinear interpolation (inside of grid) */
        *value=(1.0-a)*(1.0-b)*d[0]+a*(1.0-b)*d[1]+(1.0-a)*b*d[2]+a*b*d[3];
        *rms  =(1.0-a)*(1.0-b)*r[0]+a*(1.0-b)*r[1]+(1.0-a)*b*r[2]+a*b*r[3];
    }
    /* nearest-neighbour extrapolation (outside of grid) */
    else if (a<=0.5&&b<=0.5&&d[0]>0.0) {*value=d[0]; *rms=r[0];}
    else if (a> 0.5&&b<=0.5&&d[1]>0.0) {*value=d[1]; *rms=r[1];}
    else if (a<=0.5&&b> 0.5&&d[2]>0.0) {*value=d[2]; *rms=r[2];}
    else if (a> 0.5&&b> 0.5&&d[3]>0.0) {*value=d[3]; *rms=r[3];}
    else {
        i=0;
        for (n=0;n<4;n++) if (d[n]>0.0) {i++; *value+=d[n]; *rms+=r[n];}
        if(i==0) return 0;
        *value/=i; *rms/=i;
    }
    return 1;
}
/* ionosphere delay by tec grid data -----------------------------------------*/
static int iondelay(gtime_t time, const tec_t *tec, const double *pos,
    const double *azel, int opt, double *delay, double *var)
{
    const double fact=40.30E16/FREQ1/FREQ1; /* tecu->L1 iono (m) */
    double fs,posp[3]={0},vtec,rms,hion,rp;
    int i;

    trace(3,"iondelay: time=%s pos=%.1f %.1f azel=%.1f %.1f\n",time_str(time,0),
        pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,azel[1]*R2D);

    *delay=*var=0.0;

    for (i=0;i<tec->ndata[2];i++) { /* for a layer */

        hion=tec->hgts[0]+tec->hgts[2]*i;

        /* ionospheric pierce point position */
        fs=ionppp(pos,azel,tec->rb,hion,posp);

        if (opt&2) {
            /* modified single layer mapping function (M-SLM) ref [2] */
            rp=tec->rb/(tec->rb+hion)*sin(0.9782*(PI/2.0-azel[1]));
            fs=1.0/sqrt(1.0-rp*rp);
        }
        if (opt&1) {
            /* earth rotation correction (sun-fixed coordinate) */
            posp[1]+=2.0*PI*timediff(time,tec->time)/86400.0;
        }
        /* interpolate tec grid data */
        if (!interptec(tec,i,posp,&vtec,&rms)) return 0;

        *delay+=fact*fs*vtec;
        *var+=fact*fact*fs*fs*rms*rms;
    }
    trace(4,"iondelay: delay=%7.2f std=%6.2f\n",*delay,sqrt(*var));

    return 1;
}
/* ionosphere model by tec grid data -------------------------------------------
* compute ionospheric delay by tec grid data
* args   : gtime_t time     I   time (gpst)
*          nav_t  *nav      I   navigation data
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    opt       I   model option
*                                bit0: 0:earth-fixed,1:sun-fixed
*                                bit1: 0:single-layer,1:modified single-layer
*          double *delay    O   ionospheric delay (L1) (m)
*          double *var      O   ionospheric dealy (L1) variance (m^2)
* return : status (1:ok,0:error)
* notes  : before calling the function, read tec grid data by calling readtec()
*          return ok with delay=0 and var=VAR_NOTEC if el<MIN_EL or h<MIN_HGT
*-----------------------------------------------------------------------------*/
extern int iontec(gtime_t time, const nav_t *nav, const double *pos,
    const double *azel, int opt, double *delay, double *var)
{
    double dels[2],vars[2],a,tt;
    int i,stat[2];

    trace(3,"iontec  : time=%s pos=%.1f %.1f azel=%.1f %.1f\n",time_str(time,0),
        pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,azel[1]*R2D);

    if (azel[1]<MIN_EL||pos[2]<MIN_HGT) {
        *delay=0.0;
        *var=VAR_NOTEC;
        return 1;
    }
    for (i=0;i<nav->nt;i++) {
        if (timediff(nav->tec[i].time,time)>0.0) break;
    }
    if (i==0||i>=nav->nt) {
        trace(2,"%s: tec grid out of period\n",time_str(time,0));
        return 0;
    }
    if ((tt=timediff(nav->tec[i].time,nav->tec[i-1].time))==0.0) {
        trace(2,"tec grid time interval error\n");
        return 0;
    }
    /* ionospheric delay by tec grid data */
    stat[0]=iondelay(time,nav->tec+i-1,pos,azel,opt,dels  ,vars  );
    stat[1]=iondelay(time,nav->tec+i  ,pos,azel,opt,dels+1,vars+1);

    if (!stat[0]&&!stat[1]) {
        trace(2,"%s: tec grid out of area pos=%6.2f %7.2f azel=%6.1f %5.1f\n",
            time_str(time,0),pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,azel[1]*R2D);
        return 0;
    }
    if (stat[0]&&stat[1]) { /* linear interpolation by time */
        a=timediff(time,nav->tec[i-1].time)/tt;
        *delay=dels[0]*(1.0-a)+dels[1]*a;
        *var  =vars[0]*(1.0-a)+vars[1]*a;
    }
    else if (stat[0]) { /* nearest-neighbour extrapolation by time */
        *delay=dels[0];
        *var  =vars[0];
    }
    else {
        *delay=dels[1];
        *var  =vars[1];
    }
    trace(3,"iontec  : delay=%5.2f std=%5.2f\n",*delay,sqrt(*var));
    return 1;
}
