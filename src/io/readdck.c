
#include "io.h"

#define MAXRNXLEN   (16*MAXOBSTYPE+4)   /* max rinex record length */
#define MAXPOSHEAD  1024                /* max head line position */

/* compare precise clock -----------------------------------------------------*/
static int cmppclk(const void *p1, const void *p2)
{
    pclk_t *q1=(pclk_t*)p1,*q2=(pclk_t*)p2;
    double tt=timediff(q1->time,q2->time);
    return tt<-1E-9?-1:(tt>1E-9?1:q1->index-q2->index);
}
/* combine precise clock -----------------------------------------------------*/
static void combpclk(nav_t *nav)
{
    pclk_t *nav_pclk;
    int i,j,k,m;

    trace(3,"combpclk: nc=%d\n",nav->nc);

    if (nav->nc<=0) return;
    qsort(nav->pclk,nav->nc,sizeof(pclk_t),cmppclk);

    for (i=0,j=1;j<nav->nc;j++) {
        /* if i and j at the same epoch */
        if (fabs(timediff(nav->pclk[i].time,nav->pclk[j].time))<1E-9) {
            for (k=0;k<MAXSAT;k++) {
                if (nav->pclk[j].clk[k][0]==0.0) continue;
                for (m=0;m<3;m++) {
                    nav->pclk[i].clk[k][m]=nav->pclk[j].clk[k][m];
                    nav->pclk[i].std[k][m]=nav->pclk[j].std[k][m];
                }
            }
        }
        else if (++i<j) nav->pclk[i]=nav->pclk[j];
    }
    nav->nc=i+1;

    if (!(nav_pclk=(pclk_t*)realloc(nav->pclk,sizeof(pclk_t)*nav->nc))) {
        free(nav->pclk);nav->pclk=NULL;nav->nc=nav->ncmax=0;
        trace(1,"combpclk malloc error nc=%d\n",nav->nc);
        return;
    }
    nav->pclk=nav_pclk;
    nav->ncmax=nav->nc;

    trace(3,"combpclk: nc=%d\n",nav->nc);
}
/* read rinex clock ----------------------------------------------------------*/
static int readdckb(FILE *fp, const int mask, int index, nav_t *nav)
{
    pclk_t *nav_pclk;
    gtime_t time;
    double data[6]={0};
    int i,j,sat;
    char buff[MAXRNXLEN],satid[8]="";

    trace(3,"readdckb: index=%d\n",index);

    if (!nav) return 0;

    while (fgets(buff,sizeof(buff),fp)) {

        if (str2time(buff,8,26,&time)) {
            trace(2,"dcoupled clk invalid epoch (%34.34s)\n",buff);
            continue;
        }
        strncpy(satid,buff+3,4);

        /* only read AS (satellite clock) record */
        if (strncmp(buff,"AS",2)||!(sat=satid2no(satid))) continue;
        if (!(satsys(sat,NULL)&mask)) continue;

        for (i=0,j=40;i<6;i++,j+=20) data[i]=str2num(buff,j,19);

        if (nav->nc>=nav->ncmax) {
            nav->ncmax+=1024;
            if (!(nav_pclk=(pclk_t*)realloc(nav->pclk,sizeof(pclk_t)*(nav->ncmax)))) {
                trace(1,"readdckb malloc error: nmax=%d\n",nav->ncmax);
                free(nav->pclk);nav->pclk=NULL;nav->nc=nav->ncmax=0;
                return 0;
            }
            nav->pclk=nav_pclk;
        }

        if (nav->nc<=0||fabs(timediff(time,nav->pclk[nav->nc-1].time))>1E-9) {
            nav->nc++;
            nav->pclk[nav->nc-1].time=time;
            nav->pclk[nav->nc-1].index=index;
            for (i=0;i<MAXSAT;i++) for (j=0;j<3;j++) {
                nav->pclk[nav->nc-1].clk[i][j]=0.0;
                nav->pclk[nav->nc-1].std[i][j]=0.0f;
            }
        }
        /* nenoseconds to seconds */
        for (i=0;i<3;i++) {
            nav->pclk[nav->nc-1].clk[sat-1][i]=data[i];
            nav->pclk[nav->nc-1].std[sat-1][i]=(float)(data[i+3]);
        }
    }

    return nav->nc>0;
}
/* read clock file header ----------------------------------------------------*/
static int readdckh(FILE *fp, nav_t *nav)
{

    trace(3,"readdckh:\n");

    char buff[MAXRNXLEN],*label=buff+60;

    while (fgets(buff,MAXRNXLEN,fp)) {

        //if (strlen(buff)<=60) continue;

        //else if (strstr(label,"RINEX VERSION / TYPE")) {
        //    if ((type=*(buff+20)) != 'C') {
        //        trace(2,"unsupported clock type type=%c\n",type);
        //        return 0;
        //    }
        //    continue;
        //}
        //else if (strstr(label,"PGM / RUN BY / DATE")) continue;
        //else if (strstr(label,"COMMENT")) continue;
        if (strstr(label,"END OF HEADER")) return 1;
    }
    return 0;
}
/* read decoupled clock file ---------------------------------------------------*/
static int readdckf(const char *file, const int mask, int index, nav_t *nav)
{
    trace(3,"readdckf: file=%s index=%d\n",file,index);

    FILE* fp;
    int stat;

    if (!(fp=fopen(file,"r"))) {
        trace(2,"decoupled clock file open error (%s)\n",file);
        return 0;
    }
    /* read clock file header */
    if (!readdckh(fp,nav)) return 0;

    /* read clock file body */
    stat=readdckb(fp,mask,index,nav);

    fclose(fp);
    return stat;
}
/* read decoupled clock files --------------------------------------------------
* read decoupled clock file
* args   : char *file    I      file (wild-card * expanded)
*          int   mask    I      system mask (SYS_GPS for gps)
*          nav_t *nav    IO     navigation data    (NULL: no input)
* return : 0:error, 1:ok
*-----------------------------------------------------------------------------*/
extern int readdck(const char *file, const int mask, nav_t *nav)
{
    int i,n,index=0,stat=1;
    char *files[MAXEXFILE]={ 0 };

    trace(3,"readrnxc: file=%s mask=%d\n",file,mask);

    for (i=0;i<MAXEXFILE;i++) {
        if (!(files[i]=(char*)malloc(1024))) {
            for (i--;i>=0;i--) free(files[i]);
            return -1;
        }
    }
    /* expand wild-card */
    n=expath(file,files,MAXEXFILE);

    /* read dcoupled clock files */
    for (i=0;i<n;i++) {
        if (readdckf(files[i],mask,index++,nav))
            continue;
        stat=0;
        break;
    }
    for (i=0;i<MAXEXFILE;i++) free(files[i]);

    if (!stat) return 0;
    /* unique and combine ephemeris and precise clock */
    combpclk(nav);

    return 1;
}