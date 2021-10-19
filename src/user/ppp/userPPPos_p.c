/*------------------------------------------------------------------------------
* ppp.c : precise point positioning
*
*          Copyright (C) 2010-2018 by T.TAKASU, All rights reserved.
*
* options : -DIERS_MODEL  use IERS tide model
*           -DOUTSTAT_AMB output ambiguity parameters to solution status
*
* references :
*    [1] D.D.McCarthy, IERS Technical Note 21, IERS Conventions 1996, July 1996
*    [2] D.D.McCarthy and G.Petit, IERS Technical Note 32, IERS Conventions
*        2003, November 2003
*    [3] D.A.Vallado, Fundamentals of Astrodynamics and Applications 2nd ed,
*        Space Technology Library, 2004
*    [4] J.Kouba, A Guide to using International GNSS Service (IGS) products,
*        May 2009
*    [5] RTCM Paper, April 12, 2010, Proposed SSR Messages for SV Orbit Clock,
*        Code Biases, URA
*    [6] MacMillan et al., Atmospheric gradients and the VLBI terrestrial and
*        celestial reference frames, Geophys. Res. Let., 1997
*    [7] G.Petit and B.Luzum (eds), IERS Technical Note No. 36, IERS
*         Conventions (2010), 2010
*    [8] J.Kouba, A simplified yaw-attitude model for eclipsing GPS satellites,
*        GPS Solutions, 13:1-12, 2009
*    [9] F.Dilssner, GPS IIF-1 satellite antenna phase center and attitude
*        modeling, InsideGNSS, September, 2010
*    [10] F.Dilssner, The GLONASS-M satellite yaw-attitude model, Advances in
*        Space Research, 2010
*    [11] IGS MGEX (http://igs.org/mgex)
*-----------------------------------------------------------------------------*/

#include "ppp.h"

#define MAX_ITER        8              /* max number of iterations */


/* number of estimated states ------------------------------------------------*/
extern int userPPPnx_p(const prcopt_t *opt)
{
    return PNX(opt);
}
/* not estimate datum parameters -----------------------------------------------*/
static void rejEstAmbDatum(const rtk_t *rtk, const int *datum, int *ix)
{
    const prcopt_t *opt=&rtk->opt;
    int i,k,sat;

    if (datum[MAXSAT]) return;

    /* get datum sat */
    for (i=sat=0;i<MAXSAT;i++) if (datum[i]) {
        sat=i+1; break;
    }
    if (!sat) return;

    for (i=0;i<PNF(opt);i++) {
        k=PIB(sat,i+1,opt);
        ix[k]=0;
    }
}
/* precise point positioning -------------------------------------------------*/
extern void userPPPos_p(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    const prcopt_t *opt=&rtk->opt;
    double *xp,*Pp,*vl,*HL,*DL,*vp,*HP,*DP;
    double *rs,*dts,*vars,*varc,*azel,dr[3]={0};
    char str[32];
    const int nx=rtk->nx;
    int svh[MAXOBS],exc[MAXOBS]={0},datum[MAXSAT+1],stat=SOLQ_SINGLE;
    int i,j,nv,nc,info,*ix;

    time2str(obs[0].time,str,2);
    trace(3,"userPPPos_p : time=%s nx=%d n=%d\n",str,nx,n);

    rs=mat(6,n);dts=mat(3,n);vars=mat(1,n);varc=zeros(3,n);azel=zeros(2,n);

    for (i=0;i<MAXSAT;i++) for (j=0;j<opt->nf;j++) rtk->ssat[i].fix[j]=0;

    userUdStates_p(rtk,obs,n,nav);
    
    /* satellite positions and clocks */
    satposs(obs[0].time,obs,n,nav,rtk->opt.sateph,rs,dts,vars,svh);

    /* exclude measurements of eclipsing satellite (block IIA) */
    testeclipse(obs,n,nav,rs);

    /* earth tides correction */
    tidedisp(gpst2utc(obs[0].time),rtk->x,7,&nav->erp,opt->odisp[0],dr);

    /* mw ambiguity estimation of cnes */
    if (UCK(opt)) {
        userMWEst_cnes(rtk,obs,n,nav,rs,vars,svh,dr,exc);
        //free(rs);free(dts);free(vars);free(varc);free(azel);
        //return;
    }

    nv=PNF(opt)*2*n;nc=PNX(opt);
    ix=imat(nx,1);xp=mat(nx,1);Pp=mat(nx,nx);
    vl=mat(nv,1);HL=mat(nx,nv);DL=mat(nv,nv);
    vp=mat(nc,1);HP=mat(nx,nc);DP=mat(nc,nc);

    for (i=0;i<MAX_ITER;i++) {
        
        matcpy(xp,rtk->x,nx,1);
        matcpy(Pp,rtk->P,nx,nx);
        if (UCK(opt)) imatcpy(datum,rtk->datum,MAXSAT+1,1);
        
        /* prefit residuals */
        if (!(nv=userPPPRes_p(0,obs,n,rs,dts,vars,varc,svh,dr,exc,nav,xp,Pp,rtk,
            vl,HL,DL,vp,HP,DP,azel,ix,&nc))) {
            trace(2,"%s ppp (%d) no valid obs data\n",str,i+1);
            break;
        }

        if (UCK(opt)) {
            /* set ambiguity datum for decoupled model */
            userSetAmbDatum(rtk,obs,n,exc,azel,datum);
            /* not estimate datum parameters */
            rejEstAmbDatum(rtk,datum,ix);
        }

        /* measurement update */
        if ((info=lsqBlock(xp,Pp,ix,vl,HL,DL,vp,HP,DP,nx,nv,nc))) {
            trace(2,"%s ppp (%d) filter error info=%d\n",str,i+1,info);
            break;
        }

        /* get increment */
        for (j=0;j<nx;j++) if (ix[j]) xp[j]+=rtk->x[j];

        /* postfit residuals */
        if (userPPPRes_p(i+1,obs,n,rs,dts,vars,varc,svh,dr,exc,nav,xp,Pp,rtk,
            vl,HL,DL,vp,HP,DP,azel,ix,&nc)) {
            matcpy(rtk->x,xp,nx,1);
            matcpy(rtk->P,Pp,nx,nx);
            if (UCK(opt)) imatcpy(rtk->datum,datum,MAXSAT+1,1);
            stat=SOLQ_PPP;
            break;
        }
    }
    if (i>=MAX_ITER) {
        trace(2,"%s ppp (%d) iteration overflows\n",str,i);
    }
    if (stat==SOLQ_PPP) {
        //for (i=0;i<MAX_ITER;i++) {
        //    if (opt->ionoopt==IONOOPT_IFLC&&
        //        !ARIF_PPP(rtk,obs,n,exc,azel,xp,Pp)) break;
        //    if (opt->ionoopt==IONOOPT_EST&&
        //        !userAmbResol_uup(rtk,obs,n,nav,exc,azel,xp,Pp)) break;
        //    if (!userPPPRes_p(MAX_ITER+i,obs,n,rs,dts,vars,varc,svh,dr,exc,nav,xp,Pp,rtk,
        //        vl,HL,DL,vp,HP,DP,azel,ix,&nc)) {
        //        matcpy(xp,rtk->x,nx,1);
        //        matcpy(Pp,rtk->P,nx,nx);
        //        continue;
        //    }
        //    matcpy(rtk->xa,xp,nx,1);
        //    matcpy(rtk->Pa,Pp,nx,nx);
        //    stat=SOLQ_FIX;
        //    break;
        //}
        if (userAmbResol_cnes(rtk,obs,n,nav,exc,azel,xp,Pp)&&
        //if (userAmbResol_cnes2(rtk,obs,n,exc,azel,xp,Pp)&&
            userPPPRes_p(MAX_ITER+1,obs,n,rs,dts,vars,varc,svh,dr,exc,nav,xp,Pp,
                rtk,vl,HL,DL,vp,HP,DP,azel,ix,&nc)) {
            matcpy(rtk->xa,xp,nx,1);
            matcpy(rtk->Pa,Pp,nx,nx);
#if 0
            /* hold ambiguity */
            matcpy(rtk->x,xp,nx,1);
            matcpy(rtk->P,Pp,nx,nx);
#endif
            stat=SOLQ_FIX;
        }
       
        /* update solution status */
        userUdSolution(rtk,obs,n,stat);
        logRclk_p(rtk);
        logZtd_p(rtk);
        logIon_p(rtk,obs,n);
        logAmb_p(rtk,obs,n);
        logRes_p(rtk,obs,n);
    }
    free(rs);free(dts);free(vars);free(varc);free(azel);
    free(xp);free(Pp);free(ix);
    free(vl);free(HL);free(DL);
    free(vp);free(HP);free(DP);
}
void main_index()
{
    prcopt_t opt=prcopt_default;
    int i,f;

    FILE *fp=fopen("index_ppp.txt","w");

    opt.ionoopt=IONOOPT_IFLC;
    opt.modear=ARMODE_OFF;
    fprintf(fp,"PIC=%d\n",PIC(SYS_GPS,&opt));
    fprintf(fp,"PIT=%d\n",PIT(&opt));
    for (i=0;i<MAXSAT;i++) for (f=0;f<PNF(&opt);f++) {
        fprintf(fp,"PIB=%d,",PIB(i+1,f+1,&opt));
    }
    fprintf(fp,"\n\n");

    opt.ionoopt=IONOOPT_EST;
    opt.modear=ARMODE_CONT;
    fprintf(fp,"PIC=%d\n",PIC(SYS_GPS,&opt));
    fprintf(fp,"PIT=%d\n",PIT(&opt));
    for (i=0;i<MAXSAT;i++) {
        fprintf(fp,"PII=%d,",PII(i+1,&opt));
    }
    fprintf(fp,"\n");
    for (f=0;f<PNF(&opt);f++) fprintf(fp,"PID=%d,",PID(f+1,&opt));
    fprintf(fp,"\n");
    for (i=0;i<MAXSAT;i++) for (f=0;f<PNF(&opt);f++) {
        fprintf(fp,"PIB=%d,",PIB(i+1,f+1,&opt));
    }
    fprintf(fp,"\n");

    fclose(fp);
}