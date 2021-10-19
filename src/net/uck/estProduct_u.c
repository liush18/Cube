
#include "uck.h"

#define PNTPOS	0

#if PNTPOS
#include "pos/pos.h"
#endif

#define MAX_ITER	10
#define TRACE_UCK	0

/* update solution -----------------------------------------------------------*/
static void udSolution(net_t *uck, rcv_t *rcv, int n, int stat)
{
	const prcopt_t *opt=&uck->opt;
	int i,j,f,sat;

	uck->sol.stat=stat;
	for (i=0;i<MAXRCV;i++) rcv[i].ns=0;		/* reset number of valid satllites */
	for (i=0;i<MAXSAT;i++) uck->outc[i]++;	/* satellite outage count */

	/* update lock or unlock states of receivers and satellites */
	for (i=0;i<n&&i<MAXRCV;i++) {
		for (j=0;j<rcv[i].nobs&&j<MAXOBS;j++) {
			sat=rcv[i].obs[j].sat;

			for (f=0;f<UNF(opt);f++) {
				/* valid sat (1) or not (0) */
				if(!rcv[i].ssat[sat-1].vsat[f]) {
					rcv[i].ssat[sat-1].lock[f]=0;
					rcv[i].ssat[sat-1].outc[f]++;
					continue;
				}
				/* valid satellite of current receiver */
				rcv[i].ssat[sat-1].lock[f]++;
				rcv[i].ssat[sat-1].outc[f]=0;
				/* satllite of current receiver is locked */
				if (f==0) { uck->outc[sat-1]=0; rcv[i].ns++; }
			}
		}
		/* current receiver with valid satellite is locked */
		if (rcv[i].ns) { rcv[i].lock++; rcv[i].outc=0; }
		else { /* unlock rcv */
			rcv[i].nobs=0; rcv[i].lock=0; rcv[i].outc++;
			trace(2,"%s outlock rcv = %s\n",uck->time,rcv[i].sta.name);
			continue;
		}
	}

	/* update SNR */
	for (i=0;i<n&&i<MAXRCV;i++) {
		/* skip unlock receiver */
		if (rcv[i].outc) continue;
		for (j=0;j<rcv[i].nobs&&j<MAXOBS;j++) {
			sat=rcv[i].obs[j].sat;
			for (f=0;f<UNF(opt);f++) {
				rcv[i].ssat[sat-1].snr[f]=rcv[i].obs[j].SNR[f];
			}
		}
	}
}
/* exclude invalid rcv without valid sat -------------------------------------*/
extern int excInvalidRcv(rcv_t *rcv, int n, int ref)
{
	int i,j,sat;

	/* end program if no reference receiver */
	if (rcv[ref-1].nobs==0) {
		trace(2,"invalid datum rcv(%s)\n",rcv[ref-1].sta.name);
		return 0;
	}
	for (i=0;i<n&&i<MAXRCV;i++) {
		if (rcv[i].nobs==0) continue;
		for (j=0;j<rcv[i].nobs&&j<MAXOBS;j++) {
			sat=rcv[i].obs[j].sat;
			if (rcv[i].datum[sat-1]>OBS_REJ) break;
		}
		/* no valid obs in this receiver */
		if (j==rcv[i].nobs||j==MAXOBS) {
			trace(2,"excInvalidRcv : no valid obs rcv (%s)\n",rcv[i].sta.name);
			rcv[i].nobs=0;
		}
	}

	return 1;
}
/* information calculation of measurement equations --------------------------*/
static void calcuMeasurement_u(net_t *uck, rcv_t *rcv, int n, const nav_t *nav)
{
	const prcopt_t *opt=&uck->opt;
	obsd_t *obs;
	ssat_t *ssat;
	double rr[3],dr[3],e[3],dants[NFREQ],dantr[NFREQ];
	int i,j,sat,sys,svh[MAXOBS];

#if PNTPOS
	char msg[128]="";
	sol_t sol,sol0={{0}};
	ssat_t ssats[MAXSAT];
#endif

	/* traverse receivers */
	for (i=0;i<n&&i<MAXRCV;i++) {
		for (j=0;j<MAXSAT;j++) {
			rcv[i].ssat[j].azel[0]=rcv[i].ssat[j].azel[1]=0.0;
		}
		if (rcv[i].nobs==0) continue;
		obs=rcv[i].obs;

#if PNTPOS
		sol=sol0;
		if (!pntpos(obs,rcv[i].nobs,nav,opt,&sol,NULL,ssats,msg)) {
			trace(2,"%s spp err rcv(%s) msg(%s)\n",uck->time,rcv[i].sta.name,msg);
			rcv[i].nobs=0;
			continue;
		}
#endif
		/* satellite position */
#if 1
		satposs(obs[0].time,obs,rcv[i].nobs,nav,1,rcv[i].rs,rcv[i].dts,
			rcv[i].varrs,svh);
#else
		sat_poss(obs[0].time,rcv[i].sta.name,obs,rcv[i].nobs,nav,
			rcv[i].rs,NULL,rcv[i].varrs,NULL,svh);
#endif

		/* exclude measurements of eclipsing satellite (block IIA) DEBUG */
		//testeclipse(obs,rcv[i].nobs,nav,rcv[i].rs);

		/* earth tides correction to sinex coordinates */
		tidedisp(gpst2utc(obs[0].time),rcv[i].sta.pos,7,&nav->erp,rcv[i].odisp,
			dr);

		for (j=0;j<3;j++) rr[j]=rcv[i].sta.pos[j]+dr[j];
		ecef2pos(rr,rcv[i].pos);

#if TRACE_UCK
		if (!i) {
			trace(2,"%s rs=\n",uck->time);
			tracemat(2,rcv[i].rs,6,MAXOBS,18,3);
			trace(2,"%s rr=",uck->time);
			tracemat(2,rcv[i].sta.pos,1,3,18,5);
			trace(2,"%s dr=",uck->time);
			tracemat(2,dr,1,3,8,3);
		}
#endif

		/* every obs(sat1) for the station */
		for (j=0;j<rcv[i].nobs&&j<MAXOBS;j++) {

			sat=obs[j].sat;
			ssat=&rcv[i].ssat[sat-1];

			if (rcv[i].datum[sat-1]==OBS_REJ) continue;

#if PNTPOS
			if (!ssats[sat-1].vs) {
				rcv[i].datum[sat-1]=OBS_REJ;
				continue;
			}
#endif

			if ((rcv[i].r[j]=geodist(rcv[i].rs+j*6,rr,e))<=0.0||
				satazel(rcv[i].pos,e,ssat->azel)<opt->elmin) {
				rcv[i].datum[sat-1]=OBS_REJ;
				continue;
			}
			/* relativity correction */
			relativity(rcv[i].rs+j*6,rr,rcv[i].r+j);

			/* exclude */
			if (!(sys=satsys(sat,NULL))||
				satexclude(obs[j].sat,rcv[i].varrs[j],svh[j],opt)) {
				rcv[i].datum[sat-1]=OBS_REJ;
				continue;
			}

			/* satellite and receiver antenna model */
			satantpcv(rcv[i].rs+j*6,rr,nav->pcvs+sat-1,dants);
			antmodel(&rcv[i].pcvr,rcv[i].sta.del,ssat->azel,1,dantr);

			/* phase windup model */
			windupcorr(rcv[i].time,rcv[i].rs+j*6,rr,&ssat->phw);

			/* corrected phase and code measurement */
			if (!corrObsMeasurement(obs+j,nav,ssat->azel,opt,dantr,dants,
				ssat->phw,rcv[i].L[j],rcv[i].P[j],&rcv[i].Lc[j],
				&rcv[i].Pc[j])) {
				rcv[i].datum[sat-1]=OBS_REJ;
				continue;
			}

			rcv[i].datum[sat-1]=OBS_USE;

#if TRACE_UCK
			if (!i) {
				trace(2,"%s sat=%d r=%f strp=%f mw=%f\n",uck->time,sat,
					rcv[i].r[j],rcv[i].m_h[j]*rcv[i].zhd,rcv[i].m_w[j]);
				trace(2,"%s dants=",uck->time);
				tracemat(2,dants,1,NFREQ,8,3);
				trace(2,"%s dantr=",uck->time);
				tracemat(2,dantr,1,NFREQ,8,3);
				trace(2,"%s phw=%f, Lc=%f Pc=%f\n",uck->time,ssat->phw,
					rcv[i].Lc[j],rcv[i].Pc[j]);
			}
#endif
		}
		detCycleSlip(&uck->opt,rcv[i].ssat,obs,rcv[i].nobs,nav,rcv[i].sta.name,
			rcv[i].tt,&rcv[i].jumpc);
	}
}
/* estimate uu clock products --------------------------------------------------
* return : (0:error, 1:ok)
* notes  : flag: epoch flag (0:first,1:mid,2:last)
*			nv : the number of observation equations
*				 2 equations for each frequency of each satellite (2*UNF*nobs)
*				 1 satellite clock datum constraint				  (1)
*				 1 phase bias datum constraint for each frequency (UNF)
*			nc : the number of constraint equations
*				 random-walk for satellite phase biases			(UNSPB)
*				 random-walk for trop/iono/rev biases			((UNT+UNI+UNRPB)*n)
*				 random-walk or constants for ambiguity			(UNF*nobs)
* ----------------------------------------------------------------------------*/
extern int estProduct_u(net_t *uck, rcv_t *rcv, int n, const nav_t *nav, 
	int flag)
{
	const prcopt_t *opt=&uck->opt;
	double *xp,*Pp,*vl,*HL,*DL,*vp,*HP,*DP;
	const int nx=uck->nx,ref=uck->opt.ircv;
	const int uuar=uck->opt.ionoopt==IONOOPT_EST&&opt->modear>ARMODE_OFF;
	int i,j,nv,nc,SSC[MAXSAT],SRC[MAXRCV],stat=SOLQ_NONE,*ix,ar=0;

	trace(3,"estProduct_u: %s ref rcv (%s)\n",uck->time,rcv[ref-1].sta.name);

	uck->sol.stat=stat;
	/* reference rcv is required (first epoch if zmc, all epoch if rcv datum) */
	if (opt->datumType==DTM_REFRCV) { /* rcv datum */
		if (rcv[ref-1].nobs==0) {
			trace(2,"%s no ref rcv (%s)\n",uck->time,rcv[ref-1].sta.name);
			return 0;
		}
	}
	else { /* rcv datum (set rcv clock and phase bias to 0) */
		if (flag==FST_EPOCH&&rcv[ref-1].nobs==0) { /* first epoch */
			trace(2,"%s no ref rcv (%s)\n",uck->time,rcv[ref-1].sta.name);
			return 0;
		}
	}

	/* information calculation of measurement equations */
	calcuMeasurement_u(uck,rcv,n,nav);

	/* update states */
	udStates_u(uck,rcv,n,nav);

	/* set ambiguity transfer datum if half fixed (else not fixed solution) */
	if (uuar&&uck->fixedrec>(n/2)) udAmbDatum_u(opt,rcv,n);

	/* obs equations, satllite clock/phase bias/code bias constraints */
	nv=UNF(opt)*2*uck->nobs+1+UNF(opt)+MAXSAT;
	nc=uck->nx;
	ix=imat(nx,1);xp=mat(nx,1);Pp=mat(nx,nx);vl=mat(nv,1);vp=mat(nc,1);
	HL=mat(nx,nv);HP=mat(nx,nc);DL=mat(nv,nv);DP=mat(nc,nc);
	for (i=0;i<MAX_ITER*n;i++) {

		matcpy(xp,uck->x,nx,1);
		matcpy(Pp,uck->P,nx,nx);
		if (uuar) {
			imatcpy(SSC,uck->SSC,MAXSAT,1);
			imatcpy(SRC,uck->SRC,MAXRCV,1);
		}

		/* exclude invalid rcv */
		if (!excInvalidRcv(rcv,n,ref)) return 0;

		/* set or update datum */
		if (uuar) setDatum(uck,rcv,n,flag,SSC,SRC);

		/* design matrix and variances */
		if (!(nv=uckRes(0,uck,rcv,n,nav,xp,Pp,ix,vl,HL,DL,vp,HP,DP,&nc,flag))) {
			trace(2,"%s uck (%d) no valid obs data\n",uck->time,i+1);
			return 0;
		}
		/* lsq filter */
		if (lsqBlock(xp,Pp,ix,vl,HL,DL,vp,HP,DP,nx,nv,nc)) {
			trace(2,"lsq filter error\n");
			return 0;
		}
		/* get increment */
		for (j=0;j<nx;j++) if (ix[j]) xp[j]+=uck->x[j];
		/* post residuals */
		if (uckRes(i+1,uck,rcv,n,nav,xp,Pp,ix,vl,HL,DL,vp,HP,DP,&nc,flag)) {
			stat=SOLQ_PPP;
			matcpy(uck->x,xp,nx,1);
			matcpy(uck->P,Pp,nx,nx);
			if (opt->ionoopt==IONOOPT_EST) {
				imatcpy(uck->SSC,SSC,MAXSAT,1);
				imatcpy(uck->SRC,SRC,MAXRCV,1);
			}
			break;
		}
	}
	if (i>=MAX_ITER*n) {
		trace(2,"%s uck (%d) iteration overflows\n",uck->time,i);
		return 0;
	}
	if (stat==SOLQ_PPP) {
		if (uuar) {
			if (uckAmbResol(uck,rcv,n,xp,Pp)&&
				uckRes(MAX_ITER*n+1,uck,rcv,n,nav,xp,Pp,ix,vl,HL,DL,
					vp,HP,DP,&nc,flag)) {
				ar=1;	/* count fixed ambiguities in current epoch */
			}
			/* track but not fixed solution though ar ok */
			udTrackAmb_u(uck,rcv,n,ar);
			/* fixed solution */
			if (uck->fixedrec>(n/2)) {
				stat=SOLQ_FIX;
				if (ar) {
					matcpy(uck->x,xp,uck->nx,1);
					matcpy(uck->P,Pp,uck->nx,uck->nx);
				}
			}
			else trace(2,"%s not enough fixed rcvs num=%d\n",uck->time,
				uck->fixedrec);
		}
		udSolution(uck,rcv,n,stat);
	}

	free(ix);free(xp);free(Pp);
	free(vl);free(HL);free(DL);
	free(vp);free(HP);free(DP);
	return 1;
}