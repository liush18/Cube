/*------------------------------------------------------------------------------
* readbias.c : read carrier-phase bias for AR in un-differenced and uncombined
*-----------------------------------------------------------------------------*/

#include "io.h"

#define MAXRNXLEN   (16*MAXOBSTYPE+4)   /* max rinex record length */
#define MAXPOSHEAD  1024                /* max head line position */

/* compare carrier-phase bias ------------------------------------------------*/
static int cmpBias(const void *p1, const void *p2)
{
    bias_t *q1=(bias_t*)p1,*q2=(bias_t*)p2;
    double tt=timediff(q1->time,q2->time);
    return tt<-1E-9?-1:(tt>1E-9?1:q1->index-q2->index);
}
/* combine carrier-phase bias ------------------------------------------------*/
static void combineBias(nav_t *nav)
{
    bias_t *nav_bias;
    int i,j,k;

    trace(3,"combineBias: nb=%d\n",nav->nb);

    if (nav->nb<=0) return;
    qsort(nav->bias,nav->nb,sizeof(bias_t),cmpBias);

    for (i=0,j=1;j<nav->nb;j++) {
        if (fabs(timediff(nav->bias[i].time,nav->bias[j].time))<1E-9) {
            for (k=0;k<MAXSAT;k++) {
                if (nav->bias[j].bias[k][0]==0.0) continue;
                nav->bias[i].bias[k][0]=nav->bias[j].bias[k][0];
                nav->bias[i].std [k][0]=nav->bias[j].std [k][0];
            }
        }
        else if (++i<j) nav->bias[i]=nav->bias[j];
    }
    nav->nb=i+1;

    if (!(nav_bias=(bias_t*)realloc(nav->bias,sizeof(bias_t)*nav->nb))) {
        free(nav->bias);nav->bias=NULL;nav->nb=nav->ncmax=0;
        trace(1,"combpclk malloc error nc=%d\n",nav->nb);
        return;
    }
    nav->bias=nav_bias;
    nav->nbmax=nav->nb;

    trace(4,"combineBias: nb=%d\n",nav->nb);
}
/* read bias file body -------------------------------------------------------*/
static int readBiasBody(FILE *fp, const int mask, int index, nav_t *nav)
{
    bias_t* nav_bias;
    gtime_t time;
    double data[4],ep[6];
    int i,j,sat;
    char buff[MAXRNXLEN],satid[8]="";

    trace(3,"readBiasBody: index=%d\n",index);

    while (fgets(buff,sizeof(buff),fp)) {

        if (sscanf(buff,"%lf/%lf/%lf %lf:%lf:%lf %s",ep,ep+1,ep+2,
            ep+3,ep+4,ep+5,satid)<7) {
            trace(2,"invalid carrier-phase bias buff=%s\n",buff);
            return 0;
        }
        time=epoch2time(ep);

        /* only read satellite bias record */
        if (!(sat=satid2no(satid))) continue;

        /* skip system not in mask */
        if (!(satsys(sat,NULL)&mask)) continue;

        for (i=0,j=30;i<4;i++,j+=12) data[i]=str2num(buff,j,12);

        if (nav->nb>=nav->nbmax) {
            nav->nbmax+=1024;
            if (!(nav_bias=(bias_t*)realloc(nav->bias,sizeof(bias_t)*(nav->nbmax)))) {
                trace(1,"readBiasBody malloc error: nmax=%d\n",nav->nbmax);
                free(nav->bias);nav->bias=NULL;nav->nb=nav->nbmax=0;
                return -1;
            }
            nav->bias=nav_bias;
        }
        if (nav->nb<=0||fabs(timediff(time,nav->bias[nav->nb-1].time))>1E-9) {
            nav->nb++;
            nav->bias[nav->nb-1].time=time;
            nav->bias[nav->nb-1].index=index;
            for (i=0;i<MAXSAT;i++) for (j=0;j<NFREQ;j++) {
                nav->bias[nav->nb-1].bias[i][j]=0.0;
                nav->bias[nav->nb-1].std [i][j]=0.0f;
            }
        }
        /* only two frequencies now */
        nav->bias[nav->nb-1].bias[sat-1][0]=data[0];
        nav->bias[nav->nb-1].bias[sat-1][1]=data[1];
        nav->bias[nav->nb-1].std [sat-1][0]=(float)data[2];
        nav->bias[nav->nb-1].std [sat-1][1]=(float)data[3];
    }
    return nav->nb>0;
}
/* read single bias file -----------------------------------------------------*/
static int readBiasFile(const char *file, const int mask, int index, nav_t *nav)
{
    trace(3,"readBiasFile: file=%s index=%d\n",file,index);

    FILE* fp;
    int stat;

    if (!(fp=fopen(file,"r"))) {
        trace(2,"bias file open error: %s\n",file);
        return 0;
    }

    /* read bias file body */
    stat=readBiasBody(fp,mask,index,nav);

    fclose(fp);
    return stat;
}
/* read carrier phase bias ---------------------------------------------------
* read carrier phase bias for AR
* args   : char       *file    I      file (wild-card * expanded)
*          const int  mask     I      system mask
*          nav_t      *nav     IO     navigation data
* return : number of biases
*-----------------------------------------------------------------------------*/
extern int readBias(const char *file, const int mask, nav_t *nav)
{
    int i,n,index=0,stat=1;
    char* files[MAXEXFILE]={ 0 };

    trace(3,"readbias: file=%s mask=%d\n",file,mask);

    for (i=0;i<MAXEXFILE;i++) {
        if (!(files[i]=(char*)malloc(1024))) {
            for (i--;i>=0;i--) free(files[i]);
            return -1;
        }
    }
    /* expand wild-card */
    n=expath(file,files,MAXEXFILE);

    /* read bias files */
    for (i=0;i<n;i++) {
        if (readBiasFile(files[i],mask,index++,nav))
            continue;
        stat=0;
        break;
    }
    for (i=0;i<MAXEXFILE;i++) free(files[i]);

    if (!stat) return 0;

    /* unique and combine carrier-phase bias */
    combineBias(nav);

    return nav->nb;
}