
#include "io.h"


/* compare ephemeris ---------------------------------------------------------*/
static int cmpeph(const void *p1, const void *p2)
{
    eph_t *q1=(eph_t *)p1,*q2=(eph_t *)p2;
    return q1->ttr.time!=q2->ttr.time?(int)(q1->ttr.time-q2->ttr.time):
        (q1->toe.time!=q2->toe.time?(int)(q1->toe.time-q2->toe.time):
            q1->sat-q2->sat);
}
/* sort and unique ephemeris -------------------------------------------------*/
static void uniqeph(nav_t *nav)
{
    eph_t *nav_eph;
    int i,j;

    trace(3,"uniqeph: n=%d\n",nav->n);

    if (nav->n<=0) return;

    qsort(nav->eph,nav->n,sizeof(eph_t),cmpeph);

    for (i=1,j=0;i<nav->n;i++) {
        if (nav->eph[i].sat!=nav->eph[j].sat||
            nav->eph[i].iode!=nav->eph[j].iode) {
            nav->eph[++j]=nav->eph[i];
        }
    }
    nav->n=j+1;

    if (!(nav_eph=(eph_t *)realloc(nav->eph,sizeof(eph_t)*nav->n))) {
        trace(1,"uniqeph malloc error n=%d\n",nav->n);
        free(nav->eph); nav->eph=NULL; nav->n=nav->nmax=0;
        return;
    }
    nav->eph=nav_eph;
    nav->nmax=nav->n;

    trace(4,"uniqeph: n=%d\n",nav->n);
}
/* compare glonass ephemeris -------------------------------------------------*/
static int cmpgeph(const void *p1, const void *p2)
{
    geph_t *q1=(geph_t *)p1,*q2=(geph_t *)p2;
    return q1->tof.time!=q2->tof.time?(int)(q1->tof.time-q2->tof.time):
        (q1->toe.time!=q2->toe.time?(int)(q1->toe.time-q2->toe.time):
            q1->sat-q2->sat);
}
/* sort and unique glonass ephemeris -----------------------------------------*/
static void uniqgeph(nav_t *nav)
{
    geph_t *nav_geph;
    int i,j;

    trace(3,"uniqgeph: ng=%d\n",nav->ng);

    if (nav->ng<=0) return;

    qsort(nav->geph,nav->ng,sizeof(geph_t),cmpgeph);

    for (i=j=0;i<nav->ng;i++) {
        if (nav->geph[i].sat!=nav->geph[j].sat||
            nav->geph[i].toe.time!=nav->geph[j].toe.time||
            nav->geph[i].svh!=nav->geph[j].svh) {
            nav->geph[++j]=nav->geph[i];
        }
    }
    nav->ng=j+1;

    if (!(nav_geph=(geph_t *)realloc(nav->geph,sizeof(geph_t)*nav->ng))) {
        trace(1,"uniqgeph malloc error ng=%d\n",nav->ng);
        free(nav->geph); nav->geph=NULL; nav->ng=nav->ngmax=0;
        return;
    }
    nav->geph=nav_geph;
    nav->ngmax=nav->ng;

    trace(4,"uniqgeph: ng=%d\n",nav->ng);
}
/* compare sbas ephemeris ----------------------------------------------------*/
static int cmpseph(const void *p1, const void *p2)
{
    seph_t *q1=(seph_t *)p1,*q2=(seph_t *)p2;
    return q1->tof.time!=q2->tof.time?(int)(q1->tof.time-q2->tof.time):
        (q1->t0.time!=q2->t0.time?(int)(q1->t0.time-q2->t0.time):
            q1->sat-q2->sat);
}
/* sort and unique sbas ephemeris --------------------------------------------*/
static void uniqseph(nav_t *nav)
{
    seph_t *nav_seph;
    int i,j;

    trace(3,"uniqseph: ns=%d\n",nav->ns);

    if (nav->ns<=0) return;

    qsort(nav->seph,nav->ns,sizeof(seph_t),cmpseph);

    for (i=j=0;i<nav->ns;i++) {
        if (nav->seph[i].sat!=nav->seph[j].sat||
            nav->seph[i].t0.time!=nav->seph[j].t0.time) {
            nav->seph[++j]=nav->seph[i];
        }
    }
    nav->ns=j+1;

    if (!(nav_seph=(seph_t *)realloc(nav->seph,sizeof(seph_t)*nav->ns))) {
        trace(1,"uniqseph malloc error ns=%d\n",nav->ns);
        free(nav->seph); nav->seph=NULL; nav->ns=nav->nsmax=0;
        return;
    }
    nav->seph=nav_seph;
    nav->nsmax=nav->ns;

    trace(4,"uniqseph: ns=%d\n",nav->ns);
}
/* unique ephemerides ----------------------------------------------------------
* unique ephemerides in navigation data
* args   : nav_t *nav    IO     navigation data
* return : number of epochs
*-----------------------------------------------------------------------------*/
extern void uniqnav(nav_t *nav)
{
    trace(3,"uniqnav: neph=%d ngeph=%d nseph=%d\n",nav->n,nav->ng,nav->ns);

    /* unique ephemeris */
    uniqeph (nav);
    uniqgeph(nav);
    uniqseph(nav);
}
/* compare observation data -------------------------------------------------*/
static int cmpobs(const void *p1, const void *p2)
{
    obsd_t *q1=(obsd_t *)p1,*q2=(obsd_t *)p2;
    double tt=timediff(q1->time,q2->time);
    if (fabs(tt)>DTTOL) return tt<0?-1:1;
    if (q1->rcv!=q2->rcv) return (int)q1->rcv-(int)q2->rcv;
    return (int)q1->sat-(int)q2->sat;
}
/* sort and unique observation data --------------------------------------------
* sort and unique observation data by time, rcv, sat
* args   : obs_t *obs    IO     observation data
* return : number of epochs
*-----------------------------------------------------------------------------*/
extern int sortobs(obs_t *obs)
{
    int i,j,n;

    trace(3,"sortobs: nobs=%d\n",obs->n);

    if (obs->n<=0) return 0;

    qsort(obs->data,obs->n,sizeof(obsd_t),cmpobs);

    /* delete duplicated data */
    for (i=j=0;i<obs->n;i++) {
        if (obs->data[i].sat!=obs->data[j].sat||
            obs->data[i].rcv!=obs->data[j].rcv||
            timediff(obs->data[i].time,obs->data[j].time)!=0.0) {
            obs->data[++j]=obs->data[i];
        }
    }
    obs->n=j+1;

    for (i=n=0;i<obs->n;i=j,n++) {
        for (j=i+1;j<obs->n;j++) {
            if (timediff(obs->data[j].time,obs->data[i].time)>DTTOL) break;
        }
    }
    return n;
}
/* screen by time --------------------------------------------------------------
* screening by time start, time end, and time interval
* args   : gtime_t time  I      time
*          gtime_t ts    I      time start (ts.time==0:no screening by ts)
*          gtime_t te    I      time end   (te.time==0:no screening by te)
*          double  tint  I      time interval (s) (0.0:no screen by tint)
* return : 1:on condition, 0:not on condition
*-----------------------------------------------------------------------------*/
extern int screent(gtime_t time, gtime_t ts, gtime_t te, double tint)
{
    return (tint<=0.0||fmod(time2gpst(time,NULL)+DTTOL,tint)<=DTTOL*2.0)&&
        (ts.time==0||timediff(time,ts)>=-DTTOL)&&
        (te.time==0||timediff(time,te)<  DTTOL);
}
/* set system mask -----------------------------------------------------------*/
extern int set_sysmask(const char *opt)
{
    const char *p;
    int mask=SYS_NONE;

    if (!(p=strstr(opt,"-SYS="))) return SYS_ALL;

    for (p+=5;*p&&*p!=' ';p++) {
        switch (*p) {
        case 'G': mask|=SYS_GPS; break;
        case 'R': mask|=SYS_GLO; break;
        case 'E': mask|=SYS_GAL; break;
        case 'J': mask|=SYS_QZS; break;
        case 'C': mask|=SYS_CMP; break;
        case 'I': mask|=SYS_IRN; break;
        case 'S': mask|=SYS_SBS; break;
        }
    }
    return mask;
}
/* read RINEX file header ----------------------------------------------------*/
extern int readrnxh(FILE *fp, double *ver, char *type, int *sys, int *tsys,
    char tobs[][MAXOBSTYPE][4], nav_t *nav, sta_t *sta)
{
    double bias;
    char buff[MAXRNXLEN],*label=buff+60;
    int i=0,block=0,sat;

    trace(3,"readrnxh:\n");

    *ver=2.10; *type=' '; *sys=SYS_GPS;

    while (fgets(buff,MAXRNXLEN,fp)) {

        if (strlen(buff)<=60) {
            continue;
        }
        else if (strstr(label,"RINEX VERSION / TYPE")) {
            *ver=str2num(buff,0,9);
            *type=*(buff+20);

            /* satellite system */
            switch (*(buff+40)) {
                case ' ':
                case 'G': *sys=SYS_GPS;  *tsys=TSYS_GPS; break;
                case 'R': *sys=SYS_GLO;  *tsys=TSYS_UTC; break;
                case 'E': *sys=SYS_GAL;  *tsys=TSYS_GAL; break; /* v.2.12 */
                case 'S': *sys=SYS_SBS;  *tsys=TSYS_GPS; break;
                case 'J': *sys=SYS_QZS;  *tsys=TSYS_QZS; break; /* v.3.02 */
                case 'C': *sys=SYS_CMP;  *tsys=TSYS_CMP; break; /* v.2.12 */
                case 'I': *sys=SYS_IRN;  *tsys=TSYS_IRN; break; /* v.3.03 */
                case 'M': *sys=SYS_NONE; *tsys=TSYS_GPS; break; /* mixed */
                default :
                    trace(2,"not supported satellite system: %c\n",*(buff+40));
                    break;
            }
            continue;
        }
        else if (strstr(label,"PGM / RUN BY / DATE")) {
            continue;
        }
        else if (strstr(label,"COMMENT")) {

            /* read cnes wl satellite fractional biases */
            if (strstr(buff,"WIDELANE SATELLITE FRACTIONNAL BIASES")) block=1;
            else if (block) {
                /* cnes grg clock */
                if (!strncmp(buff,"WL",2)&&(sat=satid2no(buff+3))&&
                    sscanf(buff+40,"%lf",&bias)==1) {
                    nav->wlbias[sat-1]=bias;
                }
            }

            continue;
        }
        switch (*type) { /* file type */
            case 'O': decode_obsh(fp,buff,*ver,tsys,tobs,nav,sta); break;
            case 'N': decode_navh (buff,nav); break;
            case 'G': decode_gnavh(buff,nav); break;
            case 'H': decode_hnavh(buff,nav); break;
            case 'J': decode_navh (buff,nav); break; /* extension */
            case 'L': decode_navh (buff,nav); break; /* extension */
        }
        if (strstr(label,"END OF HEADER")) return 1;

        if (++i>=MAXPOSHEAD&&*type==' ') break; /* no RINEX file */
    }
    return 0;
}
/* read RINEX file -----------------------------------------------------------*/
static int readrnxfp(FILE *fp, gtime_t ts, gtime_t te, double tint,
    const char *opt, int flag, int index, char *type, obs_t *obs, nav_t *nav, 
    sta_t *sta)
{
    double ver;
    int sys,tsys=TSYS_GPS;
    char tobs[NUMSYS][MAXOBSTYPE][4]={{""}};

    trace(3,"readrnxfp: flag=%d index=%d\n",flag,index);

    /* read RINEX file header */
    if (!readrnxh(fp,&ver,type,&sys,&tsys,tobs,nav,sta)) return 0;

    /* flag=0:except for clock,1:clock */
    if ((!flag&&*type=='C')||(flag&&*type!='C')) return 0;

    /* read RINEX file body */
    switch (*type) {
        case 'O': return readrnxobs(fp,ts,te,tint,opt,index,ver,&tsys,tobs,obs,
            sta);
        case 'N': return readrnxnav(fp,opt,ver,sys    ,nav);
        case 'G': return readrnxnav(fp,opt,ver,SYS_GLO,nav);
        case 'H': return readrnxnav(fp,opt,ver,SYS_SBS,nav);
        case 'J': return readrnxnav(fp,opt,ver,SYS_QZS,nav); /* extension */
        case 'L': return readrnxnav(fp,opt,ver,SYS_GAL,nav); /* extension */
        case 'C': return readrnxclk(fp,opt,index,nav);
    }
    trace(2,"unsupported rinex type ver=%.2f type=%c\n",ver,*type);
    return 0;
}
/* initialize station parameter ----------------------------------------------*/
static void init_sta(sta_t *sta)
{
    int i;
    *sta->name   ='\0';
    *sta->marker ='\0';
    *sta->antdes ='\0';
    *sta->antsno ='\0';
    *sta->rectype='\0';
    *sta->recver ='\0';
    *sta->recsno ='\0';
    sta->antsetup=sta->itrf=sta->deltype=0;
    for (i=0;i<3;i++) sta->pos[i]=0.0;
    for (i=0;i<3;i++) sta->del[i]=0.0;
    sta->hgt=0.0;
}
/* uncompress and read RINEX file --------------------------------------------*/
extern int readrnxfile(const char *file, gtime_t ts, gtime_t te, double tint,
    const char *opt, int flag, int index, char *type, obs_t *obs, nav_t *nav, 
    sta_t *sta)
{
    FILE *fp;
    int cstat,stat;
    char tmpfile[1024];

    trace(3,"readrnxfile: file=%s flag=%d index=%d\n",file,flag,index);

    if (sta) init_sta(sta);

    /* uncompress file */
    if ((cstat=rtk_uncompress(file,tmpfile))<0) {
        trace(2,"rinex file uncompact error: %s\n",file);
        return 0;
    }
    if (!(fp=fopen(cstat?tmpfile:file,"r"))) {
        trace(2,"rinex file open error: %s\n",cstat?tmpfile:file);
        return 0;
    }
    /* read RINEX file */
    stat=readrnxfp(fp,ts,te,tint,opt,flag,index,type,obs,nav,sta);

    fclose(fp);

    /* delete temporary file */
    if (cstat) remove(tmpfile);

    return stat;
}
/* read RINEX OBS and NAV files ------------------------------------------------
* read RINEX OBS and NAV files
* args   : char *file    I      file (wild-card * expanded) ("": stdin)
*          int   rcv     I      receiver number for obs data
*         (gtime_t ts)   I      observation time start (ts.time==0: no limit)
*         (gtime_t te)   I      observation time end   (te.time==0: no limit)
*         (double tint)  I      observation time interval (s) (0:all)
*          char  *opt    I      RINEX options (see below,"": no option)
*          obs_t *obs    IO     observation data   (NULL: no input)
*          nav_t *nav    IO     navigation data    (NULL: no input)
*          sta_t *sta    IO     station parameters (NULL: no input)
* return : status (1:ok,0:no data,-1:error)
* notes  : read data are appended to obs and nav struct
*          before calling the function, obs and nav should be initialized.
*          observation data and navigation data are not sorted.
*          navigation data may be duplicated.
*          call sortobs() or uniqnav() to sort data or delete duplicated eph.
*
*          RINEX options (separated by spaces) :
*
*            -GLss[=shift]: select GPS signal ss (ss: RINEX 3 code, "1C","2W"...)
*            -RLss[=shift]: select GLO signal ss
*            -ELss[=shift]: select GAL signal ss
*            -JLss[=shift]: select QZS signal ss
*            -CLss[=shift]: select BDS signal ss
*            -ILss[=shift]: select IRN signal ss
*            -SLss[=shift]: select SBS signal ss
*
*                 shift: carrier phase shift to be added (cycle)
*            
*            -SYS=sys[,sys...]: select navigation systems
*                               (sys=G:GPS,R:GLO,E:GAL,J:QZS,C:BDS,I:IRN,S:SBS)
*
*-----------------------------------------------------------------------------*/
extern int readrnxt(const char *file, int rcv, gtime_t ts, gtime_t te,
    double tint, const char *opt, obs_t *obs, nav_t *nav, sta_t *sta)
{
    int i,n,stat=0;
    const char *p;
    char type=' ',*files[MAXEXFILE]={0};

    trace(3,"readrnxt: file=%s rcv=%d\n",file,rcv);

    if (!*file) {
        return readrnxfp(stdin,ts,te,tint,opt,0,1,&type,obs,nav,sta);
    }
    for (i=0;i<MAXEXFILE;i++) {
        if (!(files[i]=(char *)malloc(1024))) {
            for (i--;i>=0;i--) free(files[i]);
            return -1;
        }
    }
    /* expand wild-card */
    if ((n=expath(file,files,MAXEXFILE))<=0) {
        for (i=0;i<MAXEXFILE;i++) free(files[i]);
        return 0;
    }
    /* read rinex files */
    for (i=0;i<n&&stat>=0;i++) {
        stat=readrnxfile(files[i],ts,te,tint,opt,0,rcv,&type,obs,nav,sta);
    }
    /* if station name empty, set 4-char name from file head */
    if (type=='O'&&sta) {
        if (!(p=strrchr(file,FILEPATHSEP))) p=file-1;
        if (!*sta->name) setstr(sta->name,p+1,4);
    }
    for (i=0;i<MAXEXFILE;i++) free(files[i]);

    return stat;
}
extern int readrnx(const char *file, int rcv, const char *opt, obs_t *obs,
    nav_t *nav, sta_t *sta)
{
    gtime_t t={0};

    trace(3,"readrnx : file=%s rcv=%d\n",file,rcv);

    return readrnxt(file,rcv,t,t,0.0,opt,obs,nav,sta);
}