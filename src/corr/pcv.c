

#include "corr.h"

#define SQR(x)      ((x)*(x))


/* interpolate antenna phase center variation --------------------------------*/
static double interpvar(double ang, const double* var)
{
    double a=ang/5.0; /* ang=0-90 */
    int i=(int)a;
    if (i<0) return var[0]; else if (i>=18) return var[18];
    return var[i]*(1.0-a+i)+var[i+1]*(a-i);
}
/* receiver antenna model ------------------------------------------------------
* compute antenna offset by antenna phase center parameters
* args  :pcv_t *pcv       I   antenna phase center parameters
*          double *azel     I   azimuth/elevation for receiver {az,el} (rad)
*          int     opt      I   option (0:only offset,1:offset+pcv)
*          double *dant     O   range offsets for each frequency (m)
* return:none
* notes :current version does not support azimuth dependent terms
*-----------------------------------------------------------------------------*/
extern void antmodel(const pcv_t* pcv, const double* del, const double* azel,
    int opt, double* dant)
{
    double e[3], off[3], cosel=cos(azel[1]);
    int i, j;

    trace(4, "antmodel: azel=%6.1f %4.1f opt=%d\n", azel[0]*R2D, azel[1]*R2D, opt);

    e[0]=sin(azel[0])*cosel;
    e[1]=cos(azel[0])*cosel;
    e[2]=sin(azel[1]);

    for (i=0; i<NFREQ; i++) {
        for (j=0; j<3; j++) off[j]=pcv->off[i][j]+del[j];

        dant[i]=-dot(off, e, 3)+(opt?interpvar(90.0-azel[1]*R2D, pcv->var[i]):0.0);
    }
    trace(5, "antmodel: dant=%6.3f %6.3f\n", dant[0], dant[1]);
}
/* satellite antenna model ------------------------------------------------------
* compute satellite antenna phase center parameters
* args  :pcv_t *pcv       I   antenna phase center parameters
*          double nadir     I   nadir angle for satellite (rad)
*          double *dant     O   range offsets for each frequency (m)
* return:none
*-----------------------------------------------------------------------------*/
static void antmodel_s(const pcv_t* pcv, double nadir, double* dant)
{
    int i;

    trace(4, "antmodel_s: nadir=%6.1f\n", nadir*R2D);

    for (i=0; i<NFREQ; i++) {
        dant[i]=interpvar(nadir*R2D*5.0, pcv->var[i]);
    }
    trace(5, "antmodel_s: dant=%6.3f %6.3f\n", dant[0], dant[1]);
}
/* satellite antenna phase center variation ----------------------------------*/
extern void satantpcv(const double* rs, const double* rr, const pcv_t* pcv,
    double* dant)
{
    double ru[3], rz[3], eu[3], ez[3], nadir, cosa;
    int i;

    for (i=0; i<3; i++) {
        ru[i]=rr[i]-rs[i];
        rz[i]=-rs[i];
    }
    if (!normv3(ru, eu)||!normv3(rz, ez)) return;

    cosa=dot(eu, ez, 3);
    cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);
    nadir=acos(cosa);

    antmodel_s(pcv, nadir, dant);
}
/* search antenna parameter ----------------------------------------------------
* read satellite antenna phase center position
* args  :int    sat         I   satellite number (0: receiver antenna)
*          char   *type       I   antenna type for receiver antenna
*          gtime_t time       I   time to search parameters
*          pcvs_t *pcvs       IO  antenna parameters
* return:antenna parameter (NULL: no antenna)
*-----------------------------------------------------------------------------*/
extern pcv_t* searchpcv(int sat, const char *type, gtime_t time,
    const pcvs_t *pcvs)
{
    pcv_t* pcv;
    char buff[MAXANT],*types[2],*p;
    int i, j, n=0;

    trace(3, "searchpcv: sat=%2d type=%s\n", sat, type);

    if (sat) { /* search satellite antenna */
        for (i=0; i<pcvs->n; i++) {
            pcv=pcvs->pcv+i;
            if (pcv->sat != sat) continue;
            if (pcv->ts.time != 0&&timediff(pcv->ts, time)>0.0) continue;
            if (pcv->te.time != 0&&timediff(pcv->te, time)<0.0) continue;
            return pcv;
        }
    }
    else {
        strcpy(buff, type);
        for (p=strtok(buff, " "); p&&n<2; p=strtok(NULL, " ")) types[n++]=p;
        if (n<=0) return NULL;

        /* search receiver antenna with radome at first */
        for (i=0; i<pcvs->n; i++) {
            pcv=pcvs->pcv+i;
            for (j=0; j<n; j++) if (!strstr(pcv->type, types[j])) break;
            if (j>=n) return pcv;
        }
        /* search receiver antenna without radome */
        for (i=0; i<pcvs->n; i++) {
            pcv=pcvs->pcv+i;
            if (strstr(pcv->type, types[0]) != pcv->type) continue;

            trace(2, "pcv without radome is used type=%s\n", type);
            return pcv;
        }
    }
    return NULL;
}
/* set antenna parameters ----------------------------------------------------*/
extern void posSetpcv(gtime_t time, prcopt_t *popt, nav_t *nav, 
    const pcvs_t *pcvs, const pcvs_t *pcvr, const sta_t *sta)
{
    pcv_t *pcv,pcv0={0};
    double pos[3],del[3];
    int i,j,mode=PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_FIXED;
    char id[64];

    /* set satellite antenna parameters */
    for (i=0;i<MAXSAT;i++) {
        nav->pcvs[i]=pcv0;
        if (!(satsys(i+1,NULL)&popt->navsys)) continue;
        if (!(pcv=searchpcv(i+1,"",time,pcvs))) {
            satno2id(i+1,id);
            trace(3,"no satellite antenna pcv: %s\n",id);
            continue;
        }
        nav->pcvs[i]=*pcv;
    }
    for (i=0;i<(mode?2:1);i++) {
        popt->pcvr[i]=pcv0;
        if (!strcmp(popt->anttype[i],"*")) { /* set by station parameters */
            strcpy(popt->anttype[i],sta[i].antdes);
            if (sta[i].deltype==1) { /* xyz */
                if (norm(sta[i].pos,3)>0.0) {
                    ecef2pos(sta[i].pos,pos);
                    ecef2enu(pos,sta[i].del,del);
                    for (j=0;j<3;j++) popt->antdel[i][j]=del[j];
                }
            }
            else { /* enu */
                for (j=0;j<3;j++) popt->antdel[i][j]=sta[i].del[j];
            }
        }
        if (!(pcv=searchpcv(0,popt->anttype[i],time,pcvr))) {
            trace(2,"no receiver antenna pcv: %s\n",popt->anttype[i]);
            *popt->anttype[i]='\0';
            continue;
        }
        strcpy(popt->anttype[i],pcv->type);
        popt->pcvr[i]=*pcv;
    }
}