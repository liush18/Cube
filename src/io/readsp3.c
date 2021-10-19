#include "io.h"

/* satellite code to satellite system ----------------------------------------*/
static int code2sys(char code)
{
    if (code=='G'||code==' ') return SYS_GPS;
    if (code=='R') return SYS_GLO;
    if (code=='E') return SYS_GAL; /* SP3-d */
    if (code=='J') return SYS_QZS; /* SP3-d */
    if (code=='C') return SYS_CMP; /* SP3-d */
    if (code=='I') return SYS_IRN; /* SP3-d */
    if (code=='L') return SYS_LEO; /* SP3-d */
    return SYS_NONE;
}
/* read SP3 header -----------------------------------------------------------*/
static int readsp3h(FILE *fp, gtime_t *time, char *type, int *sats,
    double *bfact, char *tsys)
{
    int i,j,k=0,ns=0,sys,prn;
    char buff[1024];

    trace(3,"readsp3h:\n");

    for (i=0;;i++) {
        if (!fgets(buff,sizeof(buff),fp)) break;

        if (i==0) {
            *type=buff[2];
            if (str2time(buff,3,28,time)) return 0;
        }
        else if (!strncmp(buff,"+ ",2)) { /* satellite id */
            if (ns==0) {
                ns=(int)str2num(buff,4,2);
            }
            for (j=0;j<17&&k<ns;j++) {
                sys=code2sys(buff[9+3*j]);
                prn=(int)str2num(buff,10+3*j,2);
                if (k<MAXSAT) sats[k++]=satno(sys,prn);
            }
        }
        else if (!strncmp(buff,"++",2)) { /* orbit accuracy */
            continue;
        }
        else if (!strncmp(buff,"%c",2)) { /* time system */
            strncpy(tsys,buff+9,3); tsys[3]='\0';
        }
        else if (!strncmp(buff,"%f",2)&&bfact[0]==0.0) { /* fp base number */
            bfact[0]=str2num(buff, 3,10);
            bfact[1]=str2num(buff,14,12);
        }
        else if (!strncmp(buff,"%i",2)) {
            continue;
        }
        else if (!strncmp(buff,"/*",2)) { /* comment */
            continue;
        }
        else if (!strncmp(buff,"* ",2)) { /* first record */
                                          /* roll back file pointer */
            fseek(fp,-(long)strlen(buff),SEEK_CUR);
            break;
        }
    }
    return ns;
}
/* add precise ephemeris -----------------------------------------------------*/
static int addpeph(nav_t *nav, peph_t *peph)
{
    peph_t *nav_peph;

    if (nav->ne>=nav->nemax) {
        nav->nemax+=256;
        if (!(nav_peph=(peph_t *)realloc(nav->peph,sizeof(peph_t)*nav->nemax))) {
            trace(1,"readsp3b malloc error n=%d\n",nav->nemax);
            free(nav->peph); nav->peph=NULL; nav->ne=nav->nemax=0;
            return 0;
        }
        nav->peph=nav_peph;
    }
    nav->peph[nav->ne++]=*peph;
    return 1;
}
/* read SP3 body -------------------------------------------------------------*/
static void readsp3b(FILE *fp, char type, int *sats, int ns, double *bfact,
    char *tsys, int index, int opt, nav_t *nav)
{
    peph_t peph;
    gtime_t time;
    double val,std,base;
    int i,j,sat,sys,prn,n=ns*(type=='P'?1:2),pred_o,pred_c,v;
    char buff[1024];

    trace(3,"readsp3b: type=%c ns=%d index=%d opt=%d\n",type,ns,index,opt);

    while (fgets(buff,sizeof(buff),fp)) {

        if (!strncmp(buff,"EOF",3)) break;

        if (buff[0]!='*'||str2time(buff,3,28,&time)) {
            trace(2,"sp3 invalid epoch %31.31s\n",buff);
            continue;
        }
        if (!strcmp(tsys,"UTC")) time=utc2gpst(time); /* utc->gpst */
        peph.time =time;
        peph.index=index;

        for (i=0;i<MAXSAT;i++) {
            for (j=0;j<4;j++) {
                peph.pos[i][j]=0.0;
                peph.std[i][j]=0.0f;
                peph.vel[i][j]=0.0;
                peph.vst[i][j]=0.0f;
            }
            for (j=0;j<3;j++) {
                peph.cov[i][j]=0.0f;
                peph.vco[i][j]=0.0f;
            }
        }
        for (i=pred_o=pred_c=v=0;i<n&&fgets(buff,sizeof(buff),fp);i++) {

            if (strlen(buff)<4||(buff[0]!='P'&&buff[0]!='V')) continue;

            sys=buff[1]==' '?SYS_GPS:code2sys(buff[1]);
            prn=(int)str2num(buff,2,2);
            if      (sys==SYS_SBS) prn+=100;
            else if (sys==SYS_QZS) prn+=192; /* extension to sp3-c */

            if (!(sat=satno(sys,prn))) continue;

            if (buff[0]=='P') {
                pred_c=strlen(buff)>=76&&buff[75]=='P';
                pred_o=strlen(buff)>=80&&buff[79]=='P';
            }
            for (j=0;j<4;j++) {

                /* read option for predicted value */
                if (j< 3&&(opt&1)&& pred_o) continue;
                if (j< 3&&(opt&2)&&!pred_o) continue;
                if (j==3&&(opt&1)&& pred_c) continue;
                if (j==3&&(opt&2)&&!pred_c) continue;

                val=str2num(buff, 4+j*14,14);
                std=str2num(buff,61+j* 3,j<3?2:3);

                if (buff[0]=='P') { /* position */
                    if (val!=0.0&&fabs(val-999999.999999)>=1E-6) {
                        peph.pos[sat-1][j]=val*(j<3?1000.0:1E-6);
                        v=1; /* valid epoch */
                    }
                    if ((base=bfact[j<3?0:1])>0.0&&std>0.0) {
                        peph.std[sat-1][j]=(float)(pow(base,std)*(j<3?1E-3:1E-12));
                    }
                }
                else if (v) { /* velocity */
                    if (val!=0.0&&fabs(val-999999.999999)>=1E-6) {
                        peph.vel[sat-1][j]=val*(j<3?0.1:1E-10);
                    }
                    if ((base=bfact[j<3?0:1])>0.0&&std>0.0) {
                        peph.vst[sat-1][j]=(float)(pow(base,std)*(j<3?1E-7:1E-16));
                    }
                }
            }
        }
        if (v) {
            if (!addpeph(nav,&peph)) return;
        }
    }
}
/* compare precise ephemeris -------------------------------------------------*/
static int cmppeph(const void *p1, const void *p2)
{
    peph_t *q1=(peph_t *)p1,*q2=(peph_t *)p2;
    double tt=timediff(q1->time,q2->time);
    return tt<-1E-9?-1:(tt>1E-9?1:q1->index-q2->index);
}
/* combine precise ephemeris -------------------------------------------------*/
extern void combpeph(nav_t *nav, int opt)
{
    int i,j,k,m;

    trace(3,"combpeph: ne=%d\n",nav->ne);

    qsort(nav->peph,nav->ne,sizeof(peph_t),cmppeph);

    if (opt&4) return;

    for (i=0,j=1;j<nav->ne;j++) {

        if (fabs(timediff(nav->peph[i].time,nav->peph[j].time))<1E-9) {

            for (k=0;k<MAXSAT;k++) {
                if (norm(nav->peph[j].pos[k],4)<=0.0) continue;
                for (m=0;m<4;m++) nav->peph[i].pos[k][m]=nav->peph[j].pos[k][m];
                for (m=0;m<4;m++) nav->peph[i].std[k][m]=nav->peph[j].std[k][m];
                for (m=0;m<4;m++) nav->peph[i].vel[k][m]=nav->peph[j].vel[k][m];
                for (m=0;m<4;m++) nav->peph[i].vst[k][m]=nav->peph[j].vst[k][m];
            }
        }
        else if (++i<j) nav->peph[i]=nav->peph[j];
    }
    nav->ne=i+1;

    trace(4,"combpeph: ne=%d\n",nav->ne);
}
/* read sp3 precise ephemeris file ---------------------------------------------
* read sp3 precise ephemeris/clock files and set them to navigation data
* args   : char   *file       I   sp3-c precise ephemeris file
*                                 (wind-card * is expanded)
*          nav_t  *nav        IO  navigation data
*          int    opt         I   options (1: only observed + 2: only predicted +
*                                 4: not combined)
* return : none
* notes  : see ref [1]
*          precise ephemeris is appended and combined
*          nav->peph and nav->ne must by properly initialized before calling the
*          function
*          only files with extensions of .sp3, .SP3, .eph* and .EPH* are read
*-----------------------------------------------------------------------------*/
extern void readsp3(const char *file, nav_t *nav, int opt)
{
    FILE *fp;
    gtime_t time={0};
    double bfact[2]={0};
    int i,j,n,ns,sats[MAXSAT]={0};
    char *efiles[MAXEXFILE],*ext,type=' ',tsys[4]="";

    trace(3,"readpephs: file=%s\n",file);

    for (i=0;i<MAXEXFILE;i++) {
        if (!(efiles[i]=(char *)malloc(1024))) {
            for (i--;i>=0;i--) free(efiles[i]);
            return;
        }
    }
    /* expand wild card in file path */
    n=expath(file,efiles,MAXEXFILE);

    for (i=j=0;i<n;i++) {
        if (!(ext=strrchr(efiles[i],'.'))) continue;

        if (!strstr(ext,".sp3")&&!strstr(ext,".SP3")&&
            !strstr(ext,".eph")&&!strstr(ext,".EPH")) continue;
        /* "r"->"rb" fix bug of fseek in MS-DOS */
        if (!(fp=fopen(efiles[i],"rb"))) { 
            trace(2,"sp3 file open error %s\n",efiles[i]);
            continue;
        }
        /* read sp3 header */
        ns=readsp3h(fp,&time,&type,sats,bfact,tsys);

        /* read sp3 body */
        readsp3b(fp,type,sats,ns,bfact,tsys,j++,opt,nav);

        fclose(fp);
    }
    for (i=0;i<MAXEXFILE;i++) free(efiles[i]);

    /* combine precise ephemeris */
    if (nav->ne>0) combpeph(nav,opt);
}