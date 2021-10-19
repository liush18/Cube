
#include "user.h"

#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

static pcvs_t pcvss={0};        /* receiver antenna parameters */
static obs_t obss={0};          /* observation data */
static nav_t navs={0};          /* navigation data */
static sta_t stas[MAXRCV];      /* station infomation */
static int nepoch=0;            /* number of observation epochs */
static int iobsu =0;            /* current rover observation data index */
static int iobsr =0;            /* current reference observation data index */
static int revs  =0;            /* analysis direction (0:forward,1:backward) */
static int aborts=0;            /* abort status */
static sol_t *solf;             /* forward solutions */
static sol_t *solb;             /* backward solutions */
static double *rbf;             /* forward base positions */
static double *rbb;             /* backward base positions */
static int isolf=0;             /* current forward solutions index */
static int isolb=0;             /* current backward solutions index */
static char proc_rov [64]="";   /* rover for current processing */
static char proc_base[64]="";   /* base station for current processing */

extern pcvs_t *userGetPcv() { return &pcvss; }
extern obs_t  *userGetObs() { return &obss;  }
extern nav_t  *userGetNav() { return &navs;  }
extern sta_t  *userGetSta() { return stas;   }
extern int userGet_revs()   { return revs; }
extern int userGet_aborts() { return aborts; }
extern int userJudge_sol()  { return solf&&solb; } 
extern void userSet_revs(int i) { revs=i; }
extern void userSet_iobs(int i) { iobsu=iobsr=i; }
extern void userSet_isol(int i) { isolf=isolb=i; }
extern void userResetPar() { iobsu=iobsr=revs=aborts=0; }
extern void userFree_sol_rb() { free(solf); free(solb); free(rbf); free(rbb); }

extern void userMallocSol()
{
    solf=(sol_t *)malloc(sizeof(sol_t)*nepoch);
    solb=(sol_t *)malloc(sizeof(sol_t)*nepoch);
    rbf=(double *)malloc(sizeof(double)*nepoch*3);
    rbb=(double *)malloc(sizeof(double)*nepoch*3);
}
extern int userCombForw(rtk_t *rtk)
{
    int i;

    if (isolf>=nepoch) return 1;
    solf[isolf]=rtk->sol;
    for (i=0;i<3;i++) rbf[i+isolf*3]=rtk->rb[i];
    isolf++;

    return 0;
}
extern int userCombBack(rtk_t *rtk)
{
    int i;

    if (isolb>=nepoch) return 1;
    solb[isolb]=rtk->sol;
    for (i=0;i<3;i++) rbb[i+isolb*3]=rtk->rb[i];
    isolb++;

    return 0;
}
extern void userSortObs() { nepoch=sortobs(&obss); }

/* show message and check break ----------------------------------------------*/
extern int userCheckbrk(const char *format,...)
{
    va_list arg;
    char buff[1024],*p=buff;
    if (!*format) return showmsg("");
    va_start(arg,format);
    p+=vsprintf(p,format,arg);
    va_end(arg);
    if (*proc_rov&&*proc_base) sprintf(p," (%s-%s)",proc_rov,proc_base);
    else if (*proc_rov ) sprintf(p," (%s)",proc_rov );
    else if (*proc_base) sprintf(p," (%s)",proc_base);
    return showmsg(buff);
}
/* search next observation data index ----------------------------------------*/
static int nextobsf(const obs_t *obs,int *i,int rcv)
{
    double tt;
    int n;

    for (;*i<obs->n;(*i)++) if (obs->data[*i].rcv==rcv) break;
    for (n=0;*i+n<obs->n;n++) {
        tt=timediff(obs->data[*i+n].time,obs->data[*i].time);
        if (obs->data[*i+n].rcv!=rcv||tt>DTTOL) break;
    }
    return n;
}
static int nextobsb(const obs_t *obs,int *i,int rcv)
{
    double tt;
    int n;

    for (;*i>=0;(*i)--) if (obs->data[*i].rcv==rcv) break;
    for (n=0;*i-n>=0;n++) {
        tt=timediff(obs->data[*i-n].time,obs->data[*i].time);
        if (obs->data[*i-n].rcv!=rcv||tt<-DTTOL) break;
    }
    return n;
}
/* input obs data, navigation messages and sbas correction -------------------*/
extern int userInputObs(obsd_t *obs, int solq, const prcopt_t *popt)
{
    gtime_t time={0};
    int i,nu,nr,n=0;

    trace(3,"infunc  : revs=%d iobsu=%d iobsr=%d\n",revs,iobsu,iobsr);

    if (0<=iobsu&&iobsu<obss.n) {
        time=obss.data[iobsu].time;
        //settime((time=obss.data[iobsu].time));
        if (userCheckbrk("processing : %s Q=%d",time_str(time,0),solq)) {
            aborts=1; showmsg("aborted"); return -1;
        }
    }
    if (!revs) { /* input forward data */
        if ((nu=nextobsf(&obss,&iobsu,1))<=0) return -1;
        //if (popt->intpref) {
        //    for (;(nr=nextobsf(&obss,&iobsr,2))>0;iobsr+=nr)
        //        if (timediff(obss.data[iobsr].time,obss.data[iobsu].time)>-DTTOL) break;
        //}
        //else {
        for (i=iobsr;(nr=nextobsf(&obss,&i,2))>0;iobsr=i,i+=nr)
            if (timediff(obss.data[i].time,obss.data[iobsu].time)>DTTOL) break;
        //}
        nr=nextobsf(&obss,&iobsr,2);
        if (nr<=0) {
            nr=nextobsf(&obss,&iobsr,2);
        }
        for (i=0;i<nu&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsu+i];
        for (i=0;i<nr&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsr+i];
        iobsu+=nu;
    }
    else { /* input backward data */
        if ((nu=nextobsb(&obss,&iobsu,1))<=0) return -1;
        //if (popt->intpref) {
        //    for (;(nr=nextobsb(&obss,&iobsr,2))>0;iobsr-=nr)
        //        if (timediff(obss.data[iobsr].time,obss.data[iobsu].time)<DTTOL) break;
        //}
        //else {
        for (i=iobsr;(nr=nextobsb(&obss,&i,2))>0;iobsr=i,i-=nr)
            if (timediff(obss.data[i].time,obss.data[iobsu].time)<-DTTOL) break;
        //}
        nr=nextobsb(&obss,&iobsr,2);
        for (i=0;i<nu&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsu-nu+1+i];
        for (i=0;i<nr&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsr-nr+1+i];
        iobsu-=nu;
    }
    return n;
}
/* validation of combined solutions ------------------------------------------*/
static int valcomb(const sol_t *solf,const sol_t *solb)
{
    double dr[3],var[3];
    int i;
    char tstr[32];

    trace(3,"valcomb :\n");

    /* compare forward and backward solution */
    for (i=0;i<3;i++) {
        dr[i]=solf->rr[i]-solb->rr[i];
        var[i]=solf->qr[i]+solb->qr[i];
    }
    for (i=0;i<3;i++) {
        if (dr[i]*dr[i]<=16.0*var[i]) continue; /* ok if in 4-sigma */

        time2str(solf->time,tstr,2);
        trace(2,"degrade fix to float: %s dr=%.3f %.3f %.3f std=%.3f %.3f %.3f\n",
            tstr+11,dr[0],dr[1],dr[2],SQRT(var[0]),SQRT(var[1]),SQRT(var[2]));
        return 0;
    }
    return 1;
}
/* combine forward/backward solutions and output results ---------------------*/
extern void userCombres(FILE *fp,const prcopt_t *popt,const solopt_t *sopt)
{
    gtime_t time={0};
    sol_t sols={{0}},sol={{0}};
    double tt,Qf[9],Qb[9],Qs[9],rbs[3]={0},rb[3]={0},rr_f[3],rr_b[3],rr_s[3];
    int i,j,k,solstatic,pri[]={0,1,2,3,4,5,1,6};

    trace(3,"combres : isolf=%d isolb=%d\n",isolf,isolb);

    solstatic=sopt->solstatic&&
        (popt->mode==PMODE_STATIC||popt->mode==PMODE_PPP_STATIC);

    for (i=0,j=isolb-1;i<isolf&&j>=0;i++,j--) {

        if ((tt=timediff(solf[i].time,solb[j].time))<-DTTOL) {
            sols=solf[i];
            for (k=0;k<3;k++) rbs[k]=rbf[k+i*3];
            j++;
        }
        else if (tt>DTTOL) {
            sols=solb[j];
            for (k=0;k<3;k++) rbs[k]=rbb[k+j*3];
            i--;
        }
        else if (solf[i].stat<solb[j].stat) {
            sols=solf[i];
            for (k=0;k<3;k++) rbs[k]=rbf[k+i*3];
        }
        else if (solf[i].stat>solb[j].stat) {
            sols=solb[j];
            for (k=0;k<3;k++) rbs[k]=rbb[k+j*3];
        }
        else {
            sols=solf[i];
            sols.time=timeadd(sols.time,-tt/2.0);

            if ((popt->mode==PMODE_KINEMA||popt->mode==PMODE_MOVEB)&&
                sols.stat==SOLQ_FIX) {

                /* degrade fix to float if validation failed */
                if (!valcomb(solf+i,solb+j)) sols.stat=SOLQ_FLOAT;
            }
            for (k=0;k<3;k++) {
                Qf[k+k*3]=solf[i].qr[k];
                Qb[k+k*3]=solb[j].qr[k];
            }
            Qf[1]=Qf[3]=solf[i].qr[3];
            Qf[5]=Qf[7]=solf[i].qr[4];
            Qf[2]=Qf[6]=solf[i].qr[5];
            Qb[1]=Qb[3]=solb[j].qr[3];
            Qb[5]=Qb[7]=solb[j].qr[4];
            Qb[2]=Qb[6]=solb[j].qr[5];

            if (popt->mode==PMODE_MOVEB) {
                for (k=0;k<3;k++) rr_f[k]=solf[i].rr[k]-rbf[k+i*3];
                for (k=0;k<3;k++) rr_b[k]=solb[j].rr[k]-rbb[k+j*3];
                if (smoother(rr_f,Qf,rr_b,Qb,3,rr_s,Qs)) continue;
                for (k=0;k<3;k++) sols.rr[k]=rbs[k]+rr_s[k];
            }
            else {
                if (smoother(solf[i].rr,Qf,solb[j].rr,Qb,3,sols.rr,Qs)) continue;
            }
            sols.qr[0]=(float)Qs[0];
            sols.qr[1]=(float)Qs[4];
            sols.qr[2]=(float)Qs[8];
            sols.qr[3]=(float)Qs[1];
            sols.qr[4]=(float)Qs[5];
            sols.qr[5]=(float)Qs[2];

            /* smoother for velocity solution */
            //if (popt->dynamics) {
            //    for (k=0;k<3;k++) {
            //        Qf[k+k*3]=solf[i].qv[k];
            //        Qb[k+k*3]=solb[j].qv[k];
            //    }
            //    Qf[1]=Qf[3]=solf[i].qv[3];
            //    Qf[5]=Qf[7]=solf[i].qv[4];
            //    Qf[2]=Qf[6]=solf[i].qv[5];
            //    Qb[1]=Qb[3]=solb[j].qv[3];
            //    Qb[5]=Qb[7]=solb[j].qv[4];
            //    Qb[2]=Qb[6]=solb[j].qv[5];
            //    if (smoother(solf[i].rr+3,Qf,solb[j].rr+3,Qb,3,sols.rr+3,Qs)) continue;
            //    sols.qv[0]=(float)Qs[0];
            //    sols.qv[1]=(float)Qs[4];
            //    sols.qv[2]=(float)Qs[8];
            //    sols.qv[3]=(float)Qs[1];
            //    sols.qv[4]=(float)Qs[5];
            //    sols.qv[5]=(float)Qs[2];
            //}
        }
        if (!solstatic) {
            outsol(fp,&sols,rbs,sopt);
        }
        else if (time.time==0||pri[sols.stat]<=pri[sol.stat]) {
            sol=sols;
            for (k=0;k<3;k++) rb[k]=rbs[k];
            if (time.time==0||timediff(sols.time,time)<0.0) {
                time=sols.time;
            }
        }
    }
    if (solstatic&&time.time!=0.0) {
        sol.time=time;
        outsol(fp,&sol,rb,sopt);
    }
}
/* free obs and nav data -----------------------------------------------------*/
extern void userFreeObs()
{
    obs_t *obs=userGetObs();

    trace(3,"userFreeObs:\n");

    nepoch=0;
    free(obs->data); obs->data=NULL; obs->n=obs->nmax=0;
}
/* output header -------------------------------------------------------------*/
static void userOutHeader(FILE *fp, const prcopt_t *popt, const solopt_t *sopt)
{
    const char *s1[]={"GPST","UTC","JST"};
    gtime_t ts,te;
    double t1,t2;
    int i,j,w1,w2;
    char s2[32],s3[32];

    trace(3,"outheader:\n");

    if (sopt->posf==SOLF_NMEA||sopt->posf==SOLF_STAT) {
        return;
    }
    if (sopt->outhead) {
        if (!*sopt->prog) {
            fprintf(fp,"%s program   : RTKLIB ver.%s\n",COMMENTH,VER_RTKLIB);
        }
        else {
            fprintf(fp,"%s program   : %s\n",COMMENTH,sopt->prog);
        }
        for (i=0;i<obss.n;i++)    if (obss.data[i].rcv==1) break;
        for (j=obss.n-1;j>=0;j--) if (obss.data[j].rcv==1) break;
        if (j<i) {fprintf(fp,"\n%s no rover obs data\n",COMMENTH); return;}
        ts=obss.data[i].time;
        te=obss.data[j].time;
        t1=time2gpst(ts,&w1);
        t2=time2gpst(te,&w2);
        if (sopt->times>=1) ts=gpst2utc(ts);
        if (sopt->times>=1) te=gpst2utc(te);
        if (sopt->times==2) ts=timeadd(ts,9*3600.0);
        if (sopt->times==2) te=timeadd(te,9*3600.0);
        time2str(ts,s2,1);
        time2str(te,s3,1);
        fprintf(fp,"%s obs start : %s %s (week%04d %8.1fs)\n",COMMENTH,s2,s1[sopt->times],w1,t1);
        fprintf(fp,"%s obs end   : %s %s (week%04d %8.1fs)\n",COMMENTH,s3,s1[sopt->times],w2,t2);
    }
    if (sopt->outopt) {
        outprcopt(fp,popt);
    }
    if (sopt->outhead||sopt->outopt) fprintf(fp,"%s\n",COMMENTH);

    outsolhead(fp,sopt);
}
/* write header to output file -----------------------------------------------*/
extern int UserOutHead(const char *outfile, const prcopt_t *popt,
    const solopt_t *sopt)
{
    FILE *fp=stdout;

    trace(3,"outhead: outfile=%s\n",outfile);

    if (*outfile) {
        createdir(outfile);

        if (!(fp=fopen(outfile,"w"))) {
            showmsg("error : open output file %s",outfile);
            return 0;
        }
    }
    /* output header */
    userOutHeader(fp,popt,sopt);

    if (*outfile) fclose(fp);
    return 1;
}
/* initialize rtk control ------------------------------------------------------
* initialize rtk control struct
* args   : rtk_t    *rtk    IO  rtk control/result struct
*          prcopt_t *opt    I   positioning options (see rtklib.h)
* return : none
*-----------------------------------------------------------------------------*/
extern void userRTKInit(rtk_t *rtk, const prcopt_t *opt, int nx)
{
    sol_t sol0={{0}};
    ambc_t ambc0={{{0}}};
    ssat_t ssat0={0};
    ambr_t ambr={{0}};
    ambt_t amb0={{0}};
    int i;

    trace(3,"rtkinit :\n");

    rtk->sol=sol0;
    for (i=0;i<6;i++) rtk->rb[i]=0.0;
    rtk->nx=nx;
    rtk->na=nx;
    rtk->nw=MAXSAT+1;
    rtk->tt=0.0;
    rtk->x=zeros(rtk->nx,1);
    rtk->P=zeros(rtk->nx,rtk->nx);
    rtk->xa=zeros(rtk->na,1);
    rtk->Pa=zeros(rtk->na,rtk->na);
    rtk->xw=zeros(rtk->nw,1);
    rtk->Pw=zeros(rtk->nw,rtk->nw);
    rtk->nfix=rtk->neb=0;
    rtk->jumpc=0;
    for (i=0;i<MAXSAT+1;i++) rtk->datum[i]=0;
    for (i=0;i<MAXSAT;i++) {
        rtk->ambc[i]=ambc0;
        rtk->ssat[i]=ssat0;
        rtk->amb [i]=amb0;
        rtk->ambr[i]=ambr;
    }
    for (i=0;i<MAXERRMSG;i++) rtk->errbuf[i]=0;
    rtk->opt=*opt;
}
/* free rtk control ------------------------------------------------------------
* free memory for rtk control struct
* args   : rtk_t    *rtk    IO  rtk control/result struct
* return : none
*-----------------------------------------------------------------------------*/
extern void userRTKFree(rtk_t *rtk)
{
    trace(3,"rtkfree :\n");

    rtk->nx=rtk->na=0;
    free(rtk->x ); rtk->x =NULL;
    free(rtk->P ); rtk->P =NULL;
    free(rtk->xa); rtk->xa=NULL;
    free(rtk->Pa); rtk->Pa=NULL;
    free(rtk->xw); rtk->xw=NULL;
    free(rtk->Pw); rtk->Pw=NULL;
}