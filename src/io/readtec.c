/*------------------------------------------------------------------------------
* ionex.c : ionex functions
*
*          Copyright (C) 2011-2013 by T.TAKASU, All rights reserved.
*
* references:
*     [1] S.Schear, W.Gurtner and J.Feltens, IONEX: The IONosphere Map EXchange
*         Format Version 1, February 25, 1998
*     [2] S.Schaer, R.Markus, B.Gerhard and A.S.Timon, Daily Global Ionosphere
*         Maps based on GPS Carrier Phase Data Routinely producted by CODE
*         Analysis Center, Proceeding of the IGS Analysis Center Workshop, 1996
*
* version : $Revision:$ $Date:$
* history : 2011/03/29 1.0 new
*           2013/03/05 1.1 change api readtec()
*                          fix problem in case of lat>85deg or lat<-85deg
*           2014/02/22 1.2 fix problem on compiled as C++
*-----------------------------------------------------------------------------*/
#include "io.h"




/* get index -----------------------------------------------------------------*/
static int getindex(double value, const double *range)
{
    if (range[2]==0.0) return 0;
    if (range[1]>0.0&&(value<range[0]||range[1]<value)) return -1;
    if (range[1]<0.0&&(value<range[1]||range[0]<value)) return -1;
    return (int)floor((value-range[0])/range[2]+0.5);
}
/* get number of items -------------------------------------------------------*/
static int nitem(const double *range)
{
    /* the number of intervals divided,
       if range[1] is replace with range[0], the result is 1, 
       meaning range[0] belongs to the first interval;
       if the parameter is range[1], meaning range[1] belongs to the last 
       interval, and the result is also the total intervals */
    return getindex(range[1],range)+1;
}
/* data index (i:lat,j:lon,k:hgt) --------------------------------------------*/
static int dataindex(int i, int j, int k, const int *ndata)
{
    if (i<0||ndata[0]<=i||j<0||ndata[1]<=j||k<0||ndata[2]<=k) return -1;
    return i+ndata[0]*(j+ndata[1]*k);
}
/* add tec data to navigation data -------------------------------------------*/
static tec_t *addtec(const double *lats, const double *lons, const double *hgts,
    double rb, nav_t *nav)
{
    tec_t *p,*nav_tec;
    gtime_t time0={0};
    int i,n,ndata[3];

    trace(3,"addtec  :\n");

    ndata[0]=nitem(lats);   /* intervals devided by latitude */
    ndata[1]=nitem(lons);   /* intervals devided by longitude */
    ndata[2]=nitem(hgts);   /* intervals devided by height */
    if (ndata[0]<=1||ndata[1]<=1||ndata[2]<=0) return NULL;

    if (nav->nt>=nav->ntmax) {
        nav->ntmax+=256;
        if (!(nav_tec=(tec_t *)realloc(nav->tec,sizeof(tec_t)*nav->ntmax))) {
            trace(1,"readionex malloc error ntmax=%d\n",nav->ntmax);
            free(nav->tec); nav->tec=NULL; nav->nt=nav->ntmax=0;
            return NULL;
        }
        nav->tec=nav_tec;
    }
    p=nav->tec+nav->nt;
    p->time=time0;
    p->rb=rb;
    for (i=0;i<3;i++) {
        p->ndata[i]=ndata[i];
        p->lats[i]=lats[i];
        p->lons[i]=lons[i];
        p->hgts[i]=hgts[i];
    }
    /* total number of data, or grids */
    n=ndata[0]*ndata[1]*ndata[2];

    if (!(p->data=(double *)malloc(sizeof(double)*n))||
        !(p->rms =(float  *)malloc(sizeof(float )*n))) {
        return NULL;
    }
    for (i=0;i<n;i++) {
        p->data[i]=0.0;
        p->rms [i]=0.0f;
    }
    nav->nt++;
    return p;
}
/* read ionex dcb aux data ----------------------------------------------------*/
static void readionexdcb(FILE *fp, double *dcb, double *rms)
{
    int i,sat;
    char buff[1024],id[32],*label;

    trace(3,"readionexdcb:\n");

    for (i=0;i<MAXSAT;i++) dcb[i]=rms[i]=0.0;

    while (fgets(buff,sizeof(buff),fp)) {
        if (strlen(buff)<60) continue;
        label=buff+60;

        if (strstr(label,"PRN / BIAS / RMS")==label) {

            strncpy(id,buff+3,3); id[3]='\0';

            if (!(sat=satid2no(id))) {
                trace(2,"ionex invalid satellite: %s\n",id);
                continue;
            }
            dcb[sat-1]=str2num(buff, 6,10);
            rms[sat-1]=str2num(buff,16,10);
        }
        else if (strstr(label,"END OF AUX DATA")==label) break;
    }
}
/* read ionex header ---------------------------------------------------------*/
static double readionexh(FILE *fp, double *lats, double *lons, double *hgts,
    double *rb, double *nexp, double *dcb, double *rms)
{
    double ver=0.0;
    char buff[1024],*label;

    trace(3,"readionexh:\n");

    while (fgets(buff,sizeof(buff),fp)) {

        if (strlen(buff)<60) continue;
        label=buff+60;

        if (strstr(label,"IONEX VERSION / TYPE")==label) {
            if (buff[20]=='I') ver=str2num(buff,0,8);
        }
        else if (strstr(label,"BASE RADIUS")==label) {
            *rb=str2num(buff,0,8);
        }
        else if (strstr(label,"HGT1 / HGT2 / DHGT")==label) {
            hgts[0]=str2num(buff, 2,6);
            hgts[1]=str2num(buff, 8,6);
            hgts[2]=str2num(buff,14,6);
        }
        else if (strstr(label,"LAT1 / LAT2 / DLAT")==label) {
            lats[0]=str2num(buff, 2,6);
            lats[1]=str2num(buff, 8,6);
            lats[2]=str2num(buff,14,6);
        }
        else if (strstr(label,"LON1 / LON2 / DLON")==label) {
            lons[0]=str2num(buff, 2,6);
            lons[1]=str2num(buff, 8,6);
            lons[2]=str2num(buff,14,6);
        }
        else if (strstr(label,"EXPONENT")==label) {
            *nexp=str2num(buff,0,6);
        }
        else if (strstr(label,"START OF AUX DATA")==label&&
            strstr(buff,"DIFFERENTIAL CODE BIASES")) {
            readionexdcb(fp,dcb,rms);
        }
        else if (strstr(label,"END OF HEADER")==label) {
            return ver;
        }
    }
    return 0.0;
}
/* read ionex body -----------------------------------------------------------*/
static int readionexb(FILE *fp, const double *lats, const double *lons,
    const double *hgts, double rb, double nexp, nav_t *nav)
{
    tec_t *p=NULL;
    gtime_t time={0};
    double lat,lon[3],hgt,x;
    int i,j,k,n,m,index,type=0;
    char buff[1024],*label=buff+60;

    trace(3,"readionexb:\n");

    while (fgets(buff,sizeof(buff),fp)) {

        if (strlen(buff)<60) continue;

        if (strstr(label,"START OF TEC MAP")==label) {
            /* p record data in current "START OF TEC MAP", type=1 for data */
            if ((p=addtec(lats,lons,hgts,rb,nav))) type=1;
        }
        else if (strstr(label,"END OF TEC MAP")==label) {
            type=0;
            p=NULL;
        }
        else if (strstr(label,"START OF RMS MAP")==label) {
            type=2; /* type=2 for rms */
            p=NULL;
        }
        else if (strstr(label,"END OF RMS MAP")==label) {
            type=0;
            p=NULL;
        }
        else if (strstr(label,"EPOCH OF CURRENT MAP")==label) {
            if (str2time(buff,0,36,&time)) {
                trace(2,"ionex epoch invalid: %-36.36s\n",buff);
                continue;
            }
            if (type==2) {
                /* find tec identical time with that of rms */
                for (i=nav->nt-1;i>=0;i--) {
                    if (fabs(timediff(time,nav->tec[i].time))>=1.0) continue;
                    p=nav->tec+i;
                    break;
                }
            }
            else if (p) p->time=time;
        }
        else if (strstr(label,"LAT/LON1/LON2/DLON/H")==label&&p) {
            lat   =str2num(buff, 2,6);
            lon[0]=str2num(buff, 8,6);
            lon[1]=str2num(buff,14,6);
            lon[2]=str2num(buff,20,6);
            hgt   =str2num(buff,26,6);

            i=getindex(lat,p->lats);    /* grid determined by latitude, from 0 */
            k=getindex(hgt,p->hgts);    /* grid determined by height, from 0 */
            n=nitem(lon);   /* number of grid determined by longitude, from 1 */
            /* each "LAT/LON1/LON2/DLON/H" label has n tec data */
            for (m=0;m<n;m++) {
                if (m%16==0&&!fgets(buff,sizeof(buff),fp)) break;
                /* grid determined by longitude, from 1 */
                j=getindex(lon[0]+lon[2]*m,p->lons);
                if ((index=dataindex(i,j,k,p->ndata))<0) continue;

                if ((x=str2num(buff,m%16*5,5))==9999.0) continue;

                if (type==1) p->data[index]=x*pow(10.0,nexp);
                else p->rms[index]=(float)(x*pow(10.0,nexp));
            }
        }
    }
    return 1;
}
/* combine tec grid data -----------------------------------------------------*/
static void combtec(nav_t *nav)
{
    tec_t tmp;
    int i,j,n=0;

    trace(3,"combtec : nav->nt=%d\n",nav->nt);

    for (i=0;i<nav->nt-1;i++) {
        for (j=i+1;j<nav->nt;j++) {
            if (timediff(nav->tec[j].time,nav->tec[i].time)<0.0) {
                tmp=nav->tec[i];
                nav->tec[i]=nav->tec[j];
                nav->tec[j]=tmp;
            }
        }
    }
    for (i=0;i<nav->nt;i++) {
        if (i>0&&timediff(nav->tec[i].time,nav->tec[n-1].time)==0.0) {
            free(nav->tec[n-1].data);
            free(nav->tec[n-1].rms );
            nav->tec[n-1]=nav->tec[i];
            continue;
        }
        nav->tec[n++]=nav->tec[i];
    }
    nav->nt=n;

    trace(4,"combtec : nav->nt=%d\n",nav->nt);
}
/* read ionex tec grid file ----------------------------------------------------
* read ionex ionospheric tec grid file
* args   : char   *file       I   ionex tec grid file
*                                 (wind-card * is expanded)
*          nav_t  *nav        IO  navigation data
*                                 nav->nt, nav->ntmax and nav->tec are modified
*          int    opt         I   read option (1: no clear of tec data,0:clear)
* return : 0:error, 1:ok
* notes  : see ref [1]
*-----------------------------------------------------------------------------*/
extern int readtec(const char *file, nav_t *nav, int opt)
{
    FILE *fp;
    double lats[3]={0},lons[3]={0},hgts[3]={0},rb=0.0,nexp=-1.0;
    double dcb[MAXSAT]={0},rms[MAXSAT]={0};
    int i,n;
    char *efiles[MAXEXFILE];

    trace(3,"readtec : file=%s\n",file);

    /* clear of tec grid data option */
    if (!opt) {
        free(nav->tec); nav->tec=NULL; nav->nt=nav->ntmax=0;
    }
    for (i=0;i<MAXEXFILE;i++) {
        if (!(efiles[i]=(char *)malloc(1024))) {
            for (i--;i>=0;i--) free(efiles[i]);
            return 0;
        }
    }
    /* expand wild card in file path */
    n=expath(file,efiles,MAXEXFILE);

    for (i=0;i<n;i++) {
        if (!(fp=fopen(efiles[i],"r"))) {
            trace(2,"ionex file open error %s\n",efiles[i]);
            return 0;
        }
        /* read ionex header */
        if (readionexh(fp,lats,lons,hgts,&rb,&nexp,dcb,rms)<=0.0) {
            trace(2,"ionex file format error %s\n",efiles[i]);
            return 0;
        }
        /* read ionex body */
        readionexb(fp,lats,lons,hgts,rb,nexp,nav);

        fclose(fp);
    }
    for (i=0;i<MAXEXFILE;i++) free(efiles[i]);

    /* combine tec grid data */
    if (nav->nt>0) combtec(nav);

    /* P1-P2 dcb */
    for (i=0;i<MAXSAT;i++) {
        nav->cbias[i][0]=CLIGHT*dcb[i]*1E-9; /* ns->m */
    }
    return 1;
}