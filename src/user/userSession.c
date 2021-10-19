
#include "user.h"


/* open procssing session ----------------------------------------------------*/
extern int userOpenSession(const prcopt_t *popt, const filopt_t *fopt, 
    nav_t *nav, pcvs_t *pcv)
{    
    gtime_t ptime[3];
    char file[1024];
    int i;

    trace(3,"userOpenSession:\n"); 

    ptime[1]=popt->ts;
    ptime[0]=timeadd(ptime[1],-86400.0);
    ptime[2]=timeadd(ptime[1], 86400.0);

    /* read ephemeris files */
    if (popt->ts.time>0) reppath(fopt->nav,file,popt->ts);
    else strcpy(file,fopt->nav);
    if (readrnx(file,0,"",NULL,nav,NULL)<0) {
        trace(1,"readnav insufficient memory\n");
        return 0;
    }
    if (nav->n<=0&&nav->ng<=0&&nav->ns<=0) {
        trace(2,"no nav data\n");
        return 0;
    }
    /* delete duplicated ephemeris */
    uniqnav(nav);

    /* read dcoupled clock files */
    if (popt->ts.time>0) reppath(fopt->clk,file,popt->ts);
    else strcpy(file,fopt->clk);
    if (popt->dck&&readdck(file,SYS_GPS,nav)<0) {
        trace(1,"readdck insufficient memory\n");
        return 0;
    }
    /* read precise clock files */
    if (!popt->dck&&!readrnxc(file,nav)) {
        trace(1,"readclk insufficient memory\n");
        return 0;
    }
    if (nav->nc<=0) {
        trace(2,"no clk data\n");
        return 0;
    }
    /* read precise ephemeris files */
#if 1
    if (popt->ts.time>0) for (i=0;i<3;i++) {
        reppath(fopt->sp3,file,ptime[i]);
        fprintf(stdout,"read: %s\n",file);
        readsp3(file,nav,0);
    }
    else {
        strcpy(file,fopt->sp3);
        fprintf(stdout,"read: %s\n",file);
        readsp3(file,nav,0);
    }
#else
    if (popt->ts.time>0) reppath(fopt->sp3,file,popt->ts);
    else strcpy(file,fopt->sp3);
    readsp3(file,nav,0);
#endif
    if (nav->ne<=0) {
        trace(2,"no sp3 data\n");
        return 0;
    }
    /* read code bias file */
    if (popt->ts.time>0) reppath(fopt->cob,file,popt->ts);
    else strcpy(file,fopt->cob);
    if (popt->scb&&!readcob(file,1,popt->nf,nav)) {
        trace(2,"read code bias file error\n");
        return 0;
    }
    /* read carrier-phase bias */
    if (popt->ts.time>0) reppath(fopt->bis,file,popt->ts);
    else strcpy(file,fopt->bis);
    if (popt->ionoopt==IONOOPT_EST&&popt->modear>ARMODE_OFF&&
        !readBias(file,1,nav)) {
        trace(2,"read carrier-phase bias error\n");
        return 0;
    }
    /* read receiver and satellite antenna parameters */
    if (*fopt->atx&&!readpcv(fopt->atx,pcv)) {
        trace(2,"sat antenna pcv read error: %s\n",fopt->atx);
        return 0;
    }
    /* use satellite L2 offset if L5 offset does not exists */
    for (i=0;i<pcv->n;i++) {
        if (norm(pcv->pcv[i].off[2],3)>0.0) continue;
        matcpy(pcv->pcv[i].off[2],pcv->pcv[i].off[1],3,1);
        matcpy(pcv->pcv[i].var[2],pcv->pcv[i].var[1],19,1);
    }

    /* read erp data */
    if (popt->ts.time>0) reppath(fopt->erp,file,popt->ts);
    else strcpy(file,fopt->erp);
    if (file&&!readerp(file,&nav->erp)) {
        trace(2,"no erp data %s\n",file);
        return 0;
    }
    /* read dcb parameters */
    if (popt->ts.time>0) reppath(fopt->dcb,file,popt->ts);
    else strcpy(file,fopt->dcb);
    if (file&&!readdcb(file,nav,NULL)) {
        trace(2,"no dcb data %s\n",file);
        return 0;
    }

    return 1;
}
/* close procssing session ---------------------------------------------------*/
extern void userCloseSession(nav_t *nav,pcvs_t *pcv)
{
    trace(3,"closeses:\n");

    free(nav->eph ); nav->eph =NULL; nav->n =nav->nmax =0;
    free(nav->geph); nav->geph=NULL; nav->ng=nav->ngmax=0;
    free(nav->seph); nav->seph=NULL; nav->ns=nav->nsmax=0;
    free(nav->peph); nav->peph=NULL; nav->ne=nav->nemax=0;
    free(nav->pclk); nav->pclk=NULL; nav->nc=nav->ncmax=0;
    free(nav->tec) ; nav->tec =NULL; nav->nt=nav->ntmax=0;
    free(nav->bias); nav->bias=NULL; nav->nb=nav->nbmax=0;
    free(nav->ucd);  nav->ucd =NULL; nav->nd=nav->ndmax=0;
    free(nav->erp.data); nav->erp.data=NULL; nav->erp.n=nav->erp.nmax=0;

    free(pcv->pcv); pcv->pcv=NULL; pcv->n=pcv->nmax=0;

    /* close solution statistics and debug trace */
    //closeOutSolution();
    //traceclose();
}