/*------------------------------------------------------------------------------
* readcob.c : read uncalibrated code delay
*-----------------------------------------------------------------------------*/

#include "io.h"

#define MAXRNXLEN   (16*MAXOBSTYPE+4)   /* max rinex record length */
#define MAXPOSHEAD  1024                /* max head line position */

/* compare ucd ---------------------------------------------------------------*/
static int cmpUCD(const void *p1, const void *p2)
{
    bias_t *q1=(bias_t*)p1,*q2=(bias_t*)p2;
    double tt=timediff(q1->time,q2->time);
    return tt<-1E-9?-1:(tt>1E-9?1:q1->index-q2->index);
}
/* combine ucd ---------------------------------------------------------------*/
static void combineUCD(nav_t *nav)
{
    bias_t *nav_bias;
    int i,j,k;

    trace(3,"combineBias: nb=%d\n",nav->nd);

    if (nav->nd<=0) return;
    qsort(nav->ucd,nav->nd,sizeof(bias_t),cmpUCD);

    for (i=0,j=1;j<nav->nd;j++) {
        if (fabs(timediff(nav->ucd[i].time,nav->ucd[j].time))<1E-9) {
            for (k=0;k<MAXSAT;k++) {
                if (nav->ucd[j].bias[k][0]==0.0) continue;
                nav->ucd[i].bias[k][0]=nav->ucd[j].bias[k][0];
                nav->ucd[i].std [k][0]=nav->ucd[j].std [k][0];
            }
        }
        else if (++i<j) nav->ucd[i]=nav->ucd[j];
    }
    nav->nd=i+1;

    if (!(nav_bias=(bias_t*)realloc(nav->ucd,sizeof(bias_t)*nav->nd))) {
        free(nav->ucd);nav->ucd=NULL;nav->nd=nav->ncmax=0;
        trace(1,"combpclk malloc error nc=%d\n",nav->nd);
        return;
    }
    nav->ucd=nav_bias;
    nav->ndmax=nav->nd;

    trace(4,"combineUCD: nb=%d\n",nav->nd);
}
/* read ucd file body --------------------------------------------------------*/
static int readUCDBody(FILE *fp, const int mask, int index, int nf, nav_t *nav)
{
    bias_t* nav_ucd;
    gtime_t time;
    double data[4],ep[6];
    int i,j,sat;
    char buff[MAXRNXLEN],satid[8]="";

    trace(3,"readUCDBody: index=%d\n",index);

    while (fgets(buff,sizeof(buff),fp)) {

        if (sscanf(buff,"%lf/%lf/%lf %lf:%lf:%lf %s",ep,ep+1,ep+2,
            ep+3,ep+4,ep+5,satid)<7) {
            trace(2,"invalid ucd buff=%s\n",buff);
            return 0;
        }
        time=epoch2time(ep);

        /* only read satellite ucd record */
        if (!(sat=satid2no(satid))) continue;

        /* skip system not in mask */
        if (!(satsys(sat,NULL)&mask)) continue;

        for (i=0,j=30;i<2*nf;i++,j+=12) data[i]=str2num(buff,j,12);

        if (nav->nd>=nav->ndmax) {
            nav->ndmax+=1024;
            if (!(nav_ucd=(bias_t*)realloc(nav->ucd,sizeof(bias_t)*(nav->ndmax)))) {
                trace(1,"readUCDBody malloc error: nmax=%d\n",nav->ndmax);
                free(nav->ucd);nav->ucd=NULL;nav->nd=nav->ndmax=0;
                return -1;
            }
            nav->ucd=nav_ucd;
        }
        if (nav->nd<=0||fabs(timediff(time,nav->ucd[nav->nd-1].time))>1E-9) {
            nav->nd++;
            nav->ucd[nav->nd-1].time=time;
            nav->ucd[nav->nd-1].index=index;
            for (i=0;i<MAXSAT;i++) for (j=0;j<NFREQ;j++) {
                nav->ucd[nav->nd-1].bias[i][j]=0.0;
                nav->ucd[nav->nd-1].std [i][j]=0.0f;
            }
        }
        for (i=0;i<nf;i++) {
            nav->ucd[nav->nd-1].bias[sat-1][i]=data[i];
            nav->ucd[nav->nd-1].std [sat-1][i]=(float)data[i+nf];
        }
    }
    return nav->nd>0;
}
/* read single cob file ------------------------------------------------------*/
static int readcobfile(const char *file, const int mask, int index, int nf,
    nav_t *nav)
{
    trace(3,"readcobfile: file=%s index=%d\n",file,index);

    FILE* fp;
    int stat;

    if (!(fp=fopen(file,"r"))) {
        trace(2,"ucd file open error: %s\n",file);
        return 0;
    }

    /* read ucd file body */
    stat=readUCDBody(fp,mask,index,nf,nav);

    fclose(fp);
    return stat;
}
/* read code delay -------------------------------------------------------------
* read uncalibrated code delay
* args   : char       *file    I      file (wild-card * expanded)
*          const int  mask     I      system mask
*          int        nf       I      number of data frequency
*          nav_t      *nav     IO     navigation data
* return : epoch number of ucd
*-----------------------------------------------------------------------------*/
extern int readcob(const char *file, int mask, int nf, nav_t *nav)
{
    int i,n,index=0,stat=1;
    char* files[MAXEXFILE]={ 0 };

    trace(3,"readcob: file=%s mask=%d\n",file,mask);

    for (i=0;i<MAXEXFILE;i++) {
        if (!(files[i]=(char*)malloc(1024))) {
            for (i--;i>=0;i--) free(files[i]);
            return -1;
        }
    }
    /* expand wild-card */
    n=expath(file,files,MAXEXFILE);

    /* read ucd files */
    for (i=0;i<n;i++) {
        if (readcobfile(files[i],mask,index++,nf,nav))
            continue;
        stat=0;
        break;
    }
    for (i=0;i<MAXEXFILE;i++) free(files[i]);

    if (!stat) return 0;

    /* unique and combine ucd */
    combineUCD(nav);

    return nav->nd;
}