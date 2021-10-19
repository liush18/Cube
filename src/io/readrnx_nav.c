
#include "io.h"

#define SQR(x)      ((x)*(x))

static const double ura_eph[]={         /* RAa values (ref [3] 20.3.3.3.1.1) */
    2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,
    3072.0,6144.0,0.0
};

/* decode RINEX NAV header ---------------------------------------------------*/
extern void decode_navh(char *buff, nav_t *nav)
{
    int i,j;
    char *label=buff+60;

    trace(4,"decode_navh:\n");

    if      (strstr(label,"ION ALPHA"           )) { /* opt ver.2 */
        if (nav) {
            for (i=0,j=2;i<4;i++,j+=12) nav->ion_gps[i]=str2num(buff,j,12);
        }
    }
    else if (strstr(label,"ION BETA"            )) { /* opt ver.2 */
        if (nav) {
            for (i=0,j=2;i<4;i++,j+=12) nav->ion_gps[i+4]=str2num(buff,j,12);
        }
    }
    else if (strstr(label,"DELTA-UTC: A0,A1,T,W")) { /* opt ver.2 */
        if (nav) {
            for (i=0,j=3;i<2;i++,j+=19) nav->utc_gps[i]=str2num(buff,j,19);
            for (;i<4;i++,j+=9) nav->utc_gps[i]=str2num(buff,j,9);
        }
    }
    else if (strstr(label,"IONOSPHERIC CORR"    )) { /* opt ver.3 */
        if (nav) {
            if (!strncmp(buff,"GPSA",4)) {
                for (i=0,j=5;i<4;i++,j+=12) nav->ion_gps[i]=str2num(buff,j,12);
            }
            else if (!strncmp(buff,"GPSB",4)) {
                for (i=0,j=5;i<4;i++,j+=12) nav->ion_gps[i+4]=str2num(buff,j,12);
            }
            else if (!strncmp(buff,"GAL",3)) {
                for (i=0,j=5;i<4;i++,j+=12) nav->ion_gal[i]=str2num(buff,j,12);
            }
            else if (!strncmp(buff,"QZSA",4)) { /* v.3.02 */
                for (i=0,j=5;i<4;i++,j+=12) nav->ion_qzs[i]=str2num(buff,j,12);
            }
            else if (!strncmp(buff,"QZSB",4)) { /* v.3.02 */
                for (i=0,j=5;i<4;i++,j+=12) nav->ion_qzs[i+4]=str2num(buff,j,12);
            }
            else if (!strncmp(buff,"BDSA",4)) { /* v.3.02 */
                for (i=0,j=5;i<4;i++,j+=12) nav->ion_cmp[i]=str2num(buff,j,12);
            }
            else if (!strncmp(buff,"BDSB",4)) { /* v.3.02 */
                for (i=0,j=5;i<4;i++,j+=12) nav->ion_cmp[i+4]=str2num(buff,j,12);
            }
            else if (!strncmp(buff,"IRNA",4)) { /* v.3.03 */
                for (i=0,j=5;i<4;i++,j+=12) nav->ion_irn[i]=str2num(buff,j,12);
            }
            else if (!strncmp(buff,"IRNB",4)) { /* v.3.03 */
                for (i=0,j=5;i<4;i++,j+=12) nav->ion_irn[i+4]=str2num(buff,j,12);
            }
        }
    }
    else if (strstr(label,"TIME SYSTEM CORR"    )) { /* opt ver.3 */
        if (nav) {
            if (!strncmp(buff,"GPUT",4)) {
                nav->utc_gps[0]=str2num(buff, 5,17);
                nav->utc_gps[1]=str2num(buff,22,16);
                nav->utc_gps[2]=str2num(buff,38, 7);
                nav->utc_gps[3]=str2num(buff,45, 5);
            }
            else if (!strncmp(buff,"GLUT",4)) {
                nav->utc_glo[0]=-str2num(buff,5,17); /* tau_C */
            }
            else if (!strncmp(buff,"GLGP",4)) {
                nav->utc_glo[1]=str2num(buff, 5,17); /* tau_GPS */
            }
            else if (!strncmp(buff,"GAUT",4)) { /* v.3.02 */
                nav->utc_gal[0]=str2num(buff, 5,17);
                nav->utc_gal[1]=str2num(buff,22,16);
                nav->utc_gal[2]=str2num(buff,38, 7);
                nav->utc_gal[3]=str2num(buff,45, 5);
            }
            else if (!strncmp(buff,"QZUT",4)) { /* v.3.02 */
                nav->utc_qzs[0]=str2num(buff, 5,17);
                nav->utc_qzs[1]=str2num(buff,22,16);
                nav->utc_qzs[2]=str2num(buff,38, 7);
                nav->utc_qzs[3]=str2num(buff,45, 5);
            }
            else if (!strncmp(buff,"BDUT",4)) { /* v.3.02 */
                nav->utc_cmp[0]=str2num(buff, 5,17);
                nav->utc_cmp[1]=str2num(buff,22,16);
                nav->utc_cmp[2]=str2num(buff,38, 7);
                nav->utc_cmp[3]=str2num(buff,45, 5);
            }
            else if (!strncmp(buff,"SBUT",4)) { /* v.3.02 */
                nav->utc_sbs[0]=str2num(buff, 5,17);
                nav->utc_sbs[1]=str2num(buff,22,16);
                nav->utc_sbs[2]=str2num(buff,38, 7);
                nav->utc_sbs[3]=str2num(buff,45, 5);
            }
            else if (!strncmp(buff,"IRUT",4)) { /* v.3.03 */
                nav->utc_irn[0]=str2num(buff, 5,17);
                nav->utc_irn[1]=str2num(buff,22,16);
                nav->utc_irn[2]=str2num(buff,38, 7);
                nav->utc_irn[3]=str2num(buff,45, 5);
                nav->utc_irn[8]=0.0; /* A2 */
            }
        }
    }
    else if (strstr(label,"LEAP SECONDS"        )) { /* opt */
        if (nav) {
            nav->utc_gps[4]=str2num(buff, 0,6);
            nav->utc_gps[7]=str2num(buff, 6,6);
            nav->utc_gps[5]=str2num(buff,12,6);
            nav->utc_gps[6]=str2num(buff,18,6);
        }
    }
}
/* decode GNAV header --------------------------------------------------------*/
extern void decode_gnavh(char *buff, nav_t *nav)
{
    char *label=buff+60;

    trace(4,"decode_gnavh:\n");

    if      (strstr(label,"CORR TO SYTEM TIME"  )) {} /* opt */
    else if (strstr(label,"LEAP SECONDS"        )) {} /* opt */
}
/* decode GEO NAV header -----------------------------------------------------*/
extern void decode_hnavh(char *buff, nav_t *nav)
{
    char *label=buff+60;

    trace(4,"decode_hnavh:\n");

    if      (strstr(label,"CORR TO SYTEM TIME"  )) {} /* opt */
    else if (strstr(label,"D-UTC A0,A1,T,W,S,U" )) {} /* opt */
    else if (strstr(label,"LEAP SECONDS"        )) {} /* opt */
}
/* adjust time considering week handover -------------------------------------*/
static gtime_t adjweek(gtime_t t, gtime_t t0)
{
    double tt=timediff(t,t0);
    if (tt<-302400.0) return timeadd(t, 604800.0);
    if (tt> 302400.0) return timeadd(t,-604800.0);
    return t;
}
/* adjust time considering week handover -------------------------------------*/
static gtime_t adjday(gtime_t t, gtime_t t0)
{
    double tt=timediff(t,t0);
    if (tt<-43200.0) return timeadd(t, 86400.0);
    if (tt> 43200.0) return timeadd(t,-86400.0);
    return t;
}
/* URA value (m) to URA index ------------------------------------------------*/
static int uraindex(double value)
{
    int i;
    for (i=0;i<15;i++) if (ura_eph[i]>=value) break;
    return i;
}
/* Galileo SISA value (m) to SISA index --------------------------------------*/
static int sisa_index(double value)
{
    if (value<0.0 || value>6.0) return 255; /* unknown or NAPA */
    else if (value<=0.5) return (int)(value/0.01);
    else if (value<=1.0) return (int)((value-0.5)/0.02)+50;
    else if (value<=2.0) return (int)((value-1.0)/0.04)+75;
    return ((int)(value-2.0)/0.16)+100;
}
/* decode ephemeris ----------------------------------------------------------*/
static int decode_eph(double ver, int sat, gtime_t toc, const double *data,
    eph_t *eph)
{
    eph_t eph0={0};
    int sys;

    trace(4,"decode_eph: ver=%.2f sat=%2d\n",ver,sat);

    sys=satsys(sat,NULL);

    if (!(sys&(SYS_GPS|SYS_GAL|SYS_QZS|SYS_CMP|SYS_IRN))) {
        trace(3,"ephemeris error: invalid satellite sat=%2d\n",sat);
        return 0;
    }
    *eph=eph0;

    eph->sat=sat;
    eph->toc=toc;

    eph->f0=data[0];
    eph->f1=data[1];
    eph->f2=data[2];

    eph->A=SQR(data[10]); eph->e=data[ 8]; eph->i0  =data[15]; eph->OMG0=data[13];
    eph->omg =data[17]; eph->M0 =data[ 6]; eph->deln=data[ 5]; eph->OMGd=data[18];
    eph->idot=data[19]; eph->crc=data[16]; eph->crs =data[ 4]; eph->cuc =data[ 7];
    eph->cus =data[ 9]; eph->cic=data[12]; eph->cis =data[14];

    if (sys==SYS_GPS||sys==SYS_QZS) {
        eph->iode=(int)data[ 3];      /* IODE */
        eph->iodc=(int)data[26];      /* IODC */
        eph->toes=     data[11];      /* Toe (s) in GPS week */
        eph->week=(int)data[21];      /* GPS week */
        eph->toe=adjweek(gpst2time(eph->week,data[11]),toc);
        eph->ttr=adjweek(gpst2time(eph->week,data[27]),toc);

        eph->code=(int)data[20];      /* GPS: codes on L2 ch */
        eph->svh =(int)data[24];      /* SV health */
        eph->sva=uraindex(data[23]);  /* URA index (m->index) */
        eph->flag=(int)data[22];      /* GPS: L2 P data flag */

        eph->tgd[0]=   data[25];      /* TGD */
        if (sys==SYS_GPS) {
            eph->fit=data[28];        /* fit interval (h) */
        }
        else {
            eph->fit=data[28]==0.0?1.0:2.0; /* fit interval (0:1h,1:>2h) */
        }
    }
    else if (sys==SYS_GAL) { /* GAL ver.3 */
        eph->iode=(int)data[ 3];      /* IODnav */
        eph->toes=     data[11];      /* Toe (s) in Galileo week */
        eph->week=(int)data[21];      /* Galileo week = GPS week */
        eph->toe=adjweek(gpst2time(eph->week,data[11]),toc);
        eph->ttr=adjweek(gpst2time(eph->week,data[27]),toc);

        eph->code=(int)data[20];      /* data sources */
                                      /* bit 0 set: I/NAV E1-B */
                                      /* bit 1 set: F/NAV E5a-I */
                                      /* bit 2 set: F/NAV E5b-I */
                                      /* bit 8 set: af0-af2 toc are for E5a.E1 */
                                      /* bit 9 set: af0-af2 toc are for E5b.E1 */
        eph->svh =(int)data[24];      /* sv health */
                                      /* bit     0: E1B DVS */
                                      /* bit   1-2: E1B HS */
                                      /* bit     3: E5a DVS */
                                      /* bit   4-5: E5a HS */
                                      /* bit     6: E5b DVS */
                                      /* bit   7-8: E5b HS */
        eph->sva =sisa_index(data[23]); /* sisa (m->index) */

        eph->tgd[0]=   data[25];      /* BGD E5a/E1 */
        eph->tgd[1]=   data[26];      /* BGD E5b/E1 */
    }
    else if (sys==SYS_CMP) { /* BeiDou v.3.02 */
        eph->toc=bdt2gpst(eph->toc);  /* bdt -> gpst */
        eph->iode=(int)data[ 3];      /* AODE */
        eph->iodc=(int)data[28];      /* AODC */
        eph->toes=     data[11];      /* Toe (s) in BDT week */
        eph->week=(int)data[21];      /* bdt week */
        eph->toe=bdt2gpst(bdt2time(eph->week,data[11])); /* BDT -> GPST */
        eph->ttr=bdt2gpst(bdt2time(eph->week,data[27])); /* BDT -> GPST */
        eph->toe=adjweek(eph->toe,toc);
        eph->ttr=adjweek(eph->ttr,toc);

        eph->svh =(int)data[24];      /* satH1 */
        eph->sva=uraindex(data[23]);  /* URA index (m->index) */

        eph->tgd[0]=   data[25];      /* TGD1 B1/B3 */
        eph->tgd[1]=   data[26];      /* TGD2 B2/B3 */
    }
    else if (sys==SYS_IRN) { /* IRNSS v.3.03 */
        eph->iode=(int)data[ 3];      /* IODEC */
        eph->toes=     data[11];      /* Toe (s) in IRNSS week */
        eph->week=(int)data[21];      /* IRNSS week */
        eph->toe=adjweek(gpst2time(eph->week,data[11]),toc);
        eph->ttr=adjweek(gpst2time(eph->week,data[27]),toc);
        eph->svh =(int)data[24];      /* SV health */
        eph->sva=uraindex(data[23]);  /* URA index (m->index) */
        eph->tgd[0]=   data[25];      /* TGD */
    }
    if (eph->iode<0||1023<eph->iode) {
        trace(2,"rinex nav invalid: sat=%2d iode=%d\n",sat,eph->iode);
    }
    if (eph->iodc<0||1023<eph->iodc) {
        trace(2,"rinex nav invalid: sat=%2d iodc=%d\n",sat,eph->iodc);
    }
    return 1;
}
/* decode GLONASS ephemeris --------------------------------------------------*/
static int decode_geph(double ver, int sat, gtime_t toc, double *data,
    geph_t *geph)
{
    geph_t geph0={0};
    gtime_t tof;
    double tow,tod;
    int week,dow;

    trace(4,"decode_geph: ver=%.2f sat=%2d\n",ver,sat);

    if (satsys(sat,NULL)!=SYS_GLO) {
        trace(3,"glonass ephemeris error: invalid satellite sat=%2d\n",sat);
        return 0;
    }
    *geph=geph0;

    geph->sat=sat;

    /* Toc rounded by 15 min in utc */
    tow=time2gpst(toc,&week);
    toc=gpst2time(week,floor((tow+450.0)/900.0)*900);
    dow=(int)floor(tow/86400.0);

    /* time of frame in UTC */
    tod=ver<=2.99?data[2]:fmod(data[2],86400.0); /* Tod (v.2), Tow (v.3) in UTC */
    tof=gpst2time(week,tod+dow*86400.0);
    tof=adjday(tof,toc);

    geph->toe=utc2gpst(toc);   /* Toc (GPST) */
    geph->tof=utc2gpst(tof);   /* Tof (GPST) */

                               /* IODE = Tb (7bit), Tb =index of UTC+3H within current day */
    geph->iode=(int)(fmod(tow+10800.0,86400.0)/900.0+0.5);

    geph->taun=-data[0];       /* -taun */
    geph->gamn= data[1];       /* +gamman */

    geph->pos[0]=data[3]*1E3; geph->pos[1]=data[7]*1E3; geph->pos[2]=data[11]*1E3;
    geph->vel[0]=data[4]*1E3; geph->vel[1]=data[8]*1E3; geph->vel[2]=data[12]*1E3;
    geph->acc[0]=data[5]*1E3; geph->acc[1]=data[9]*1E3; geph->acc[2]=data[13]*1E3;

    geph->svh=(int)data[ 6];
    geph->frq=(int)data[10];
#if 0 /*  output dtaun instead of age */
    geph->dtaun=data[14];
#else
    geph->age=(int)data[14];
#endif    
    /* some receiver output >128 for minus frequency number */
    if (geph->frq>128) geph->frq-=256;

    if (geph->frq<MINFREQ_GLO||MAXFREQ_GLO<geph->frq) {
        trace(2,"rinex gnav invalid freq: sat=%2d fn=%d\n",sat,geph->frq);
    }
    return 1;
}
/* decode GEO ephemeris ------------------------------------------------------*/
static int decode_seph(double ver, int sat, gtime_t toc, double *data,
    seph_t *seph)
{
    seph_t seph0={0};
    int week;

    trace(4,"decode_seph: ver=%.2f sat=%2d\n",ver,sat);

    if (satsys(sat,NULL)!=SYS_SBS) {
        trace(3,"geo ephemeris error: invalid satellite sat=%2d\n",sat);
        return 0;
    }
    *seph=seph0;

    seph->sat=sat;
    seph->t0 =toc;

    time2gpst(toc,&week);
    seph->tof=adjweek(gpst2time(week,data[2]),toc);

    seph->af0=data[0];
    seph->af1=data[1];

    seph->pos[0]=data[3]*1E3; seph->pos[1]=data[7]*1E3; seph->pos[2]=data[11]*1E3;
    seph->vel[0]=data[4]*1E3; seph->vel[1]=data[8]*1E3; seph->vel[2]=data[12]*1E3;
    seph->acc[0]=data[5]*1E3; seph->acc[1]=data[9]*1E3; seph->acc[2]=data[13]*1E3;

    seph->svh=(int)data[6];
    seph->sva=uraindex(data[10]);

    return 1;
}
/* read RINEX navigation data body -------------------------------------------*/
static int readrnxnavb(FILE *fp, const char *opt, double ver, int sys,
    int *type, eph_t *eph, geph_t *geph, seph_t *seph)
{
    gtime_t toc;
    double data[64];
    int i=0,j,prn,sat=0,sp=3,mask;
    char buff[MAXRNXLEN],id[8]="",*p;

    trace(4,"readrnxnavb: ver=%.2f sys=%d\n",ver,sys);

    /* set system mask */
    mask=set_sysmask(opt);

    while (fgets(buff,MAXRNXLEN,fp)) {

        if (i==0) {

            /* decode satellite field */
            if (ver>=3.0||sys==SYS_GAL||sys==SYS_QZS) { /* ver.3 or GAL/QZS */
                sprintf(id,"%.3s",buff);
                sat=satid2no(id);
                sp=4;
                if (ver>=3.0) {
                    sys=satsys(sat,NULL);
                    if (!sys) {
                        sys=(id[0]=='S')?SYS_SBS:((id[0]=='R')?SYS_GLO:SYS_GPS);
                    }
                }
            }
            else {
                prn=(int)str2num(buff,0,2);

                if (sys==SYS_SBS) {
                    sat=satno(SYS_SBS,prn+100);
                }
                else if (sys==SYS_GLO) {
                    sat=satno(SYS_GLO,prn);
                }
                else if (93<=prn&&prn<=97) { /* extension */
                    sat=satno(SYS_QZS,prn+100);
                }
                else sat=satno(SYS_GPS,prn);
            }
            /* decode Toc field */
            if (str2time(buff+sp,0,19,&toc)) {
                trace(2,"rinex nav toc error: %23.23s\n",buff);
                return 0;
            }
            /* decode data fields */
            for (j=0,p=buff+sp+19;j<3;j++,p+=19) {
                data[i++]=str2num(p,0,19);
            }
        }
        else {
            /* decode data fields */
            for (j=0,p=buff+sp;j<4;j++,p+=19) {
                data[i++]=str2num(p,0,19);
            }
            /* decode ephemeris */
            if (sys==SYS_GLO&&i>=15) {
                if (!(mask&sys)) return 0;
                *type=1;
                return decode_geph(ver,sat,toc,data,geph);
            }
            else if (sys==SYS_SBS&&i>=15) {
                if (!(mask&sys)) return 0;
                *type=2;
                return decode_seph(ver,sat,toc,data,seph);
            }
            else if (i>=31) {
                if (!(mask&sys)) return 0;
                *type=0;
                return decode_eph(ver,sat,toc,data,eph);
            }
        }
    }
    return -1;
}
/* add ephemeris to navigation data ------------------------------------------*/
static int add_eph(nav_t *nav, const eph_t *eph)
{
    eph_t *nav_eph;

    if (nav->nmax<=nav->n) {
        nav->nmax+=1024;
        if (!(nav_eph=(eph_t *)realloc(nav->eph,sizeof(eph_t)*nav->nmax))) {
            trace(1,"decode_eph malloc error: n=%d\n",nav->nmax);
            free(nav->eph); nav->eph=NULL; nav->n=nav->nmax=0;
            return 0;
        }
        nav->eph=nav_eph;
    }
    nav->eph[nav->n++]=*eph;
    return 1;
}
static int add_geph(nav_t *nav, const geph_t *geph)
{
    geph_t *nav_geph;

    if (nav->ngmax<=nav->ng) {
        nav->ngmax+=1024;
        if (!(nav_geph=(geph_t *)realloc(nav->geph,sizeof(geph_t)*nav->ngmax))) {
            trace(1,"decode_geph malloc error: n=%d\n",nav->ngmax);
            free(nav->geph); nav->geph=NULL; nav->ng=nav->ngmax=0;
            return 0;
        }
        nav->geph=nav_geph;
    }
    nav->geph[nav->ng++]=*geph;
    return 1;
}
static int add_seph(nav_t *nav, const seph_t *seph)
{
    seph_t *nav_seph;

    if (nav->nsmax<=nav->ns) {
        nav->nsmax+=1024;
        if (!(nav_seph=(seph_t *)realloc(nav->seph,sizeof(seph_t)*nav->nsmax))) {
            trace(1,"decode_seph malloc error: n=%d\n",nav->nsmax);
            free(nav->seph); nav->seph=NULL; nav->ns=nav->nsmax=0;
            return 0;
        }
        nav->seph=nav_seph;
    }
    nav->seph[nav->ns++]=*seph;
    return 1;
}
/* read RINEX navigation data ------------------------------------------------*/
extern int readrnxnav(FILE *fp, const char *opt, double ver, int sys, nav_t *nav)
{
    eph_t eph;
    geph_t geph;
    seph_t seph;
    int stat,type;

    trace(3,"readrnxnav: ver=%.2f sys=%d\n",ver,sys);

    if (!nav) return 0;

    /* read RINEX navigation data body */
    while ((stat=readrnxnavb(fp,opt,ver,sys,&type,&eph,&geph,&seph))>=0) {

        /* add ephemeris to navigation data */
        if (stat) {
            switch (type) {
            case 1 : stat=add_geph(nav,&geph); break;
            case 2 : stat=add_seph(nav,&seph); break;
            default: stat=add_eph (nav,&eph ); break;
            }
            if (!stat) return 0;
        }
    }
    return nav->n>0||nav->ng>0||nav->ns>0;
}
