
#include "dckp.h"

#define MAXPRCDAYS  100          /* max days of continuous processing */

/* process positioning ---------------------------------------------------------
mode: 0 for forward or backward, 1 for combined 
revs: 0 for forward, 1 for backwatd 
------------------------------------------------------------------------------*/
static void procpos(FILE *fp, const prcopt_t *popt, const solopt_t *sopt, 
    int mode)
{
    gtime_t time={0};
    sol_t sol={{0}};
    rtk_t rtk;
    obsd_t obs[MAXOBS*2]; /* for rover and base */
    nav_t *nav=userGetNav();
    sta_t *sta=userGetSta();
    double rb[3]={0};
    int i,nobs,n;

    trace(3,"procpos : mode=%d\n",mode);

    userRTKInit(&rtk,popt,userPPPnx_dp(popt));
    while ((nobs=userInputObs(obs,rtk.sol.stat,popt))>=0) {

        /* exclude satellites */
        for (i=n=0;i<nobs;i++) {
            if ((satsys(obs[i].sat,NULL)&popt->navsys)&&
                popt->exsats[obs[i].sat-1]!=1) obs[n++]=obs[i];
        }
        if (n<=0) continue;

        /* get current epoch time string */
        time2str(obs[0].time,rtk.cTime,2);

        detCycleSlip(popt,rtk.ssat,obs,n,nav,sta[0].name,rtk.tt,&rtk.jumpc);
        if (!userRTKPos_dp(&rtk,obs,n,nav)) continue;

        if (mode==0) { /* forward/backward */
            outsol(fp,&rtk.sol,rtk.rb,sopt);
        }
        else if (!userGet_revs()) { /* combined-forward */
            if (userCombForw(&rtk)) return;
        }
        else { /* combined-backward */
            if (userCombBack(&rtk)) return;
        }
    }
    userRTKFree(&rtk);
}
/* open output file for append -----------------------------------------------*/
static FILE *openfile(const char *outfile)
{
    trace(3,"openfile: outfile=%s\n",outfile);
    return !*outfile?stdout:fopen(outfile,"a");
}
/* execute processing session ------------------------------------------------*/
static int execses(gtime_t ts, gtime_t te, double ti, const char *obsf,
    const char *outf, const prcopt_t *popt, const solopt_t *sopt, 
    const filopt_t *fopt, int flag)
{
    FILE *fp;
    prcopt_t popt_=*popt;
    pcvs_t *pcvs=userGetPcv();
    obs_t *obs=userGetObs();
    nav_t *nav=userGetNav();
    sta_t *sta=userGetSta();

    trace(3,"execses :\n");

    /* read obs data */
    if (readrnxt(obsf,1,ts,te,ti,popt->rnxopt[0],obs,nav,sta)<0) {
        trace(1,"readobs error\n");
        return 1;
    }
    if (obs->n<=0) {
        trace(2,"no obs data\n");
        return 1;
    }
    /* sort observation data */
    userSortObs();

    /* output solution header */
    if (!flag&&!UserOutHead(outf,popt,sopt)) {
        userFreeObs();
        return 1;
    }
    /* set antenna paramters */
    if (popt_.mode!=PMODE_SINGLE) {
        posSetpcv(obs->n>0?obs->data[0].time:timeget(),&popt_,nav,pcvs,pcvs,sta);
    }
    /* read ocean tide loading parameters */
    if (popt_.mode>PMODE_SINGLE&&*fopt->blq) {
        readotl_pos(&popt_,fopt->blq,sta);
    }

    userResetPar();
    if (popt_.mode==PMODE_SINGLE||popt_.soltype==0) {
        if ((fp=openfile(outf))) {
            procpos(fp,&popt_,sopt,0); /* forward */
            fclose(fp);
        }
    }
    else if (popt_.soltype==1) {
        if ((fp=openfile(outf))) {
            userSet_revs(1); userSet_iobs(obs->n-1);
            procpos(fp,&popt_,sopt,0); /* backward */
            fclose(fp);
        }
    }
    else { /* combined */
        userMallocSol();

        if (userJudge_sol()) {
            userSet_isol(0);
            procpos(NULL,&popt_,sopt,1); /* forward */
            userSet_revs(1); userSet_iobs(obs->n-1);
            procpos(NULL,&popt_,sopt,1); /* backward */

                                         /* combine forward/backward solutions */
            if (!userGet_aborts()&&(fp=openfile(outf))) {
                userCombres(fp,&popt_,sopt);
                fclose(fp);
            }
        }
        else showmsg("error : memory allocation");
        userFree_sol_rb();
    }
    /* free obs and nav data */
    userFreeObs();

    return userGet_aborts()?1:0;
}
/* execute processing time ---------------------------------------------------*/
static int execses_t(gtime_t ts, gtime_t te, double ti, double tu,
    const char *obsf, const char *outf,const prcopt_t *popt, 
    const solopt_t *sopt,const filopt_t *fopt)
{
    trace(3,"execses_t : ti=%.0f tu=%.0f\n",ti,tu);

    gtime_t tts,tte;
    double tunit,tss;
    int i,stat=0,week;

    if (ts.time!=0&&te.time!=0&&tu>=0.0) {
        if (timediff(te,ts)<0.0) {
            showmsg("error : no period");
            return 1;
        }
        if (tu==0.0||tu>86400.0*MAXPRCDAYS) tu=86400.0*MAXPRCDAYS;
        tunit=tu<86400.0?tu:86400.0;
        tss=tunit*(int)floor(time2gpst(ts,&week)/tunit);

        for (i=0;;i++) { /* for each periods */
            tts=gpst2time(week,tss+i*tu);
            tte=timeadd(tts,tu-DTTOL);
            if (timediff(tts,te)>0.0) break;
            if (timediff(tts,ts)<0.0) tts=ts;
            if (timediff(tte,te)>0.0) tte=te;

            if (userCheckbrk("reading    : %s",time_str(tts,0))) {
                stat=1;
                break;
            }

            /* execute processing session */
            stat=execses(tts,tte,ti,obsf,outf,popt,sopt,fopt,i);

            if (stat==1) break;
        }
    }
    else {
        /* execute processing session */
        stat=execses(ts,te,ti,obsf,outf,popt,sopt,fopt,0);
    }

    return stat;
}
/* post-processing positioning -------------------------------------------------
* post-processing positioning
* args   : gtime_t ts       I   processing start time (ts.time==0: no limit)
*        : gtime_t te       I   processing end time   (te.time==0: no limit)
*          double ti        I   processing interval  (s) (0:all)
*          double tu        I   processing unit time (s) (0:all)
*          prcopt_t *popt   I   processing options
*          solopt_t *sopt   I   solution options
*          filopt_t *fopt   I   file options
* return : status (0:ok,0>:error,1:aborted)
*
*          inputs files can include wild-cards (*). if an file includes
*          wild-cards, the wild-card expanded multiple files are used.
*       
*          ts is required if time keyword in file path need to be replaced
*
*-----------------------------------------------------------------------------*/
extern int userPostPos_dp(gtime_t ts, gtime_t te, double ti, double tu,
    const prcopt_t *popt, const solopt_t *sopt, const filopt_t *fopt)
{
    pcvs_t *pcvs=userGetPcv();
    nav_t *nav=userGetNav();
    int i,n,stat=0;
    char *files[MAXEXFILE]={0},*ext,soltype[16];
    char obsfile[1024],outfile[1024],tracefile[1024],statfile[1024];

    for (i=0;i<MAXEXFILE;i++) {
        if (!(files[i]=(char*)malloc(1024))) {
            for (i--;i>=0;i--) free(files[i]);
            return -1;
        }
    }
    /* get obs files */
    if (ts.time>0) reppath(fopt->obs,obsfile,ts);
    else strcpy(obsfile,fopt->obs);
    if (strstr(fopt->obs,".list")) n=getOFilePath(obsfile,files);
    else n=expath(obsfile,files,MAXEXFILE);
    if (n<=0) {
        for (i=0;i<MAXEXFILE;i++) free(files[i]);
        return 1;
    }
    /* name solution folder */
    if (popt->modear>ARMODE_OFF) sprintf(soltype,"_fixed");
    else sprintf(soltype,"_float");
    /* processing files one by one */
    for (i=0;i<n;i++) {
        showmsg("NO.%-3d processing obs file (%s) ...\n",i+1,files[i]);

        /* configure output path */
        strcpy(outfile,files[i]);
        if (ext=strrchr(outfile,FILEPATHSEP)) 
            *ext='\0';
        if (ext=strrchr(outfile,FILEPATHSEP)) 
            sprintf(ext+1,"%s%s","ppp_solution_dc",soltype);
        if (ext=strrchr(files[i],FILEPATHSEP))
            strcat(outfile,ext);
        if (ext=strrchr(outfile,'.')) sprintf(ext,"%s",".pos");

        /* open debug trace */
        if (sopt->trace>0) {
            strcpy(tracefile,outfile);
            strcat(tracefile,".trace");
            traceclose();
            traceopen(tracefile);
            tracelevel(sopt->trace);
        }
        /* open solution statistics */
        if (sopt->sstat>0) {
            /* open receiver clock record file */
            strcpy(statfile,outfile);
            strcat(statfile,".rclk");
            logRclkClose_dp();
            logRclkOpen_dp(statfile,popt);
            /* open ztd record file */
            strcpy(statfile,outfile);
            strcat(statfile,".ztd");
            logZtdClose_dp();
            logZtdOpen_dp(statfile,popt);
            /* open ambiguity record file */
            strcpy(statfile,outfile);
            strcat(statfile,".amb");
            logAmbClose_dp();
            logAmbOpen_dp(statfile,popt);
            /* open residual record file */
            strcpy(statfile,outfile);
            strcat(statfile,".res");
            logResClose_dp();
            logResOpen_dp(statfile,popt);
        }

        /* only read corrections in the first processing */
        if (i==0) {
            if (!userOpenSession(popt,fopt,nav,pcvs)) {
                stat=1;
                break;
            }
        }
        stat=execses_t(ts,te,ti,tu,files[i],outfile,popt,sopt,fopt);
        if (stat) trace(2,"processing file=%s error\n",files[i]);
    }

    for (i=0;i<MAXEXFILE;i++) free(files[i]);

    /* close processing session */
    userCloseSession(nav,pcvs);
    logRclkClose_dp();logResClose_dp();logZtdClose_dp();logAmbClose_dp();
    traceclose();

    return stat;
}