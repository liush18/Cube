
#include "dck.h"

/* update solution -----------------------------------------------------------*/
static void udSolution_d(net_t *dck, rcv_t *rcv, int n, int stat)
{
	const prcopt_t *opt=&dck->opt;
	const int nx=dck->nx;
	int i,j,f,sat,out[MAXSAT];

	dck->sol.stat=stat;
	for (i=0;i<MAXRCV;i++) rcv[i].ns=0;		/* reset number of valid satllites */
	for (i=0;i<MAXSAT;i++) out[i]=1;

	for (i=0;i<n&&i<MAXRCV;i++) {
		for (j=0;j<rcv[i].nobs&&j<MAXOBS;j++) {
			sat=rcv[i].obs[j].sat;

			/* valid sat or not */
			if(!rcv[i].ssat[sat-1].vsat[0]) {
				rcv[i].ssat[sat-1].lock[0]=0;
				rcv[i].ssat[sat-1].outc[0]++;
				continue;
			}
			/* valid satellite of current receiver */
			rcv[i].ssat[sat-1].lock[0]++;
			rcv[i].ssat[sat-1].outc[0]=0;
			rcv[i].ns++;
			out[sat-1]=0;
		}
		if (rcv[i].ns) { 
			rcv[i].lock++; rcv[i].outc=0; 
		}
		else { /* unlock rcv */
			netRejRcvAllSat(&rcv[i]); rcv[i].lock=0; rcv[i].outc++;
			trace(2,"%s outlock rcv=%s outc=%d\n",dck->time,rcv[i].sta.name,
				rcv[i].outc);
			continue;
		}
	}

	/* outlock satellite */
	for (j=0;j<MAXSAT;j++) {
		if (out[j]) {
			dck->outc[j]++;
			trace(2,"udSolution_d: sat=%02d outc=%d\n",j+1,dck->outc[j]);
		}
		else dck->outc[j]=0;
	}

	/* update SNR */
	for (i=0;i<n&&i<MAXRCV;i++) {
		/* skip out lock rcv */
		if (rcv[i].outc) continue;
		for (j=0;j<rcv[i].nobs&&j<MAXOBS;j++) {
			sat=rcv[i].obs[j].sat;
			for (f=0;f<DNF(opt);f++) {
				rcv[i].ssat[sat-1].snr[f]=rcv[i].obs[j].SNR[f];
			}
		}
	}
}
/* exclude invalid rcv without valid sat -------------------------------------*/
extern int excInvalidRcv_d(rcv_t *rcv, int n, int ref)
{
	int i,j,sat;

	for (i=0;i<n&&i<MAXRCV;i++) {
		if (rcv[i].nobs==0) continue;
		for (j=0;j<rcv[i].nobs&&j<MAXOBS;j++) {
			sat=rcv[i].obs[j].sat;
			if (rcv[i].datum[sat-1]>OBS_REJ) break;
		}
		if (j==rcv[i].nobs||j==MAXOBS) rcv[i].nobs=0;
	}
	if (rcv[ref-1].nobs==0) {
		trace(2,"invalid ref rcv(%s) nobs=%d\n",rcv[ref-1].sta.name,rcv[i].nobs);
		return 0;
	}
	return 1;
}
/* information calculation of measurement equations --------------------------*/
static void calcuMeasurement_d(net_t *dck, rcv_t *rcv, int n, const nav_t *nav)
{
	const prcopt_t *opt=&dck->opt;
	obsd_t *obs;
	ssat_t *ssat;
	double rr[3],dr[3],e[3],dants[NFREQ],dantr[NFREQ];
	int i,j,sat,sys,svh[MAXOBS];

	/* traverse stations */
	for (i=0;i<n&&i<MAXRCV;i++) {

		for (j=0;j<MAXSAT;j++) {
			rcv[i].ssat[j].azel[0]=rcv[i].ssat[j].azel[1]=0.0;
		}
		if (rcv[i].nobs==0) continue;
		obs=rcv[i].obs;

		/* satellite position */
		sat_poss(obs[0].time,rcv[i].sta.name,obs,rcv[i].nobs,nav,
			rcv[i].rs,NULL,rcv[i].varrs,NULL,svh);

		/* exclude measurements of eclipsing satellite (block IIA) DEBUG */
		//testeclipse(obs,rcv[i].nobs,nav,rcv[i].rs);

		/* earth tides correction to sinex coordinates */
		tidedisp(gpst2utc(obs[0].time),rcv[i].sta.pos,7,&nav->erp,rcv[i].odisp,dr);

		for (j=0;j<3;j++) rr[j]=rcv[i].sta.pos[j]+dr[j];
		ecef2pos(rr,rcv[i].pos);

		/* every obs(sat1) for the station */
		for (j=0;j<rcv[i].nobs&&j<MAXOBS;j++) {

			sat=obs[j].sat;
			ssat=&rcv[i].ssat[sat-1];
			
			if (rcv[i].datum[sat-1]==OBS_REJ) continue;

			if ((rcv[i].r[j]=geodist(rcv[i].rs+j*6,rr,e))<=0.0||
				satazel(rcv[i].pos,e,ssat->azel)<opt->elmin) {
				rcv[i].datum[sat-1]=OBS_REJ;
				continue;
			}
			/* relativity correction */
			relativity(rcv[i].rs+j*6,rr,rcv[i].r+j);

			if (!(sys=satsys(sat,NULL))||
				satexclude(obs[j].sat,rcv[i].varrs[j],svh[j],opt)) {
				rcv[i].datum[sat-1]=OBS_REJ;
				continue;
			}

			/* satellite and receiver antenna model */
			satantpcv(rcv[i].rs+j*6,rr,nav->pcvs+sat-1,dants);
			antmodel(&rcv[i].pcvr,rcv[i].sta.del,ssat->azel,1,dantr);

			/* phase windup model */
			windupcorr(rcv[i].time,rcv[i].rs+j*6,rr,&rcv[i].ssat[sat-1].phw);

			/* corrected phase and code measurement */
			if (!corrObsMeasurement(obs+j,nav,ssat->azel,opt,dantr,dants,
				ssat->phw,rcv[i].L[j],rcv[i].P[j],&rcv[i].Lc[j],&rcv[i].Pc[j])) {
				rcv[i].datum[sat-1]=OBS_REJ;
				continue;
			}
			/* mw widelane combination */
			combWideLane(rcv[i].L[j],rcv[i].P[j],rcv[i].A4+j);

			if (rcv[i].Lc[j]==0.0||rcv[i].Pc[j]==0.0||rcv[i].A4[j]==0.0) {
				rcv[i].datum[sat-1]=OBS_REJ;
				continue;
			}

			rcv[i].datum[sat-1]=OBS_USE;
		}
		detCycleSlip(opt,rcv[i].ssat,obs,rcv[i].nobs,nav,rcv[i].sta.name,
			rcv[i].tt,&rcv[i].jumpc);
#if 0
		trace(2,"varr=\n");tracemat(2,varrs,    1,rcv[i].nobs,15,3);
		trace(2,"r=\n");   tracemat(2,rcv[i].r, 1,rcv[i].nobs,15,3);
		trace(2,"L3=\n");  tracemat(2,rcv[i].L3,1,rcv[i].nobs,15,3);
		trace(2,"P3=\n");  tracemat(2,rcv[i].P3,1,rcv[i].nobs,15,3);
		trace(2,"A4=\n");  tracemat(2,rcv[i].A4,1,rcv[i].nobs,15,3);
		trace(2,"\n");
#endif
	}
}
/* estimate decoupled clock products -----------------------------------------*/
extern int estProduct_d(net_t *dck, rcv_t *rcv, int n, const nav_t *nav, int flag)
{
	const prcopt_t *opt=&dck->opt;
	double *xp,*Pp,*vl,*HL,*DL,*vp,*HP,*DP;
	const int nx=dck->nx,ref=dck->opt.ircv;
	int outr[MAXRCV]={0},outs[MAXSAT]={0},SSC[MAXSAT],SRC[MAXRCV];
	int i,j,nv,nc,tnv,stat=SOLQ_NONE,states=1,*ix,ar=0;

	trace(3,"estProduct_d: %s ref rcv (%s)\n",dck->time,rcv[ref-1].sta.name);

	dck->sol.stat=stat;
	/* reference rcv is required (first epoch if zmc, all epoch if rcv datum) */
	if (opt->datumType==DTM_REFRCV) { /* rcv datum */
		if (rcv[ref-1].nobs==0) {
			trace(2,"%s no ref rcv (%s)\n",dck->time,rcv[ref-1].sta.name);
			return 0;
		}
	}
	else { /* rcv datum (set rcv clock to 0) */
		if (flag==FST_EPOCH&&rcv[ref-1].nobs==0) { /* first epoch */
			trace(2,"%s no ref rcv (%s)\n",dck->time,rcv[ref-1].sta.name);
			return 0;
		}
	}
	
	/* information calculation of measurement equations */
	calcuMeasurement_d(dck,rcv,n,nav);

	/* update states (trop and bias) */
	udStates_d(dck,rcv,n);

	/* set ambiguity transfer datum if half fixed (else not fixed solution) */
	if (dck->fixedrec>(n/2)) udAmbDatum_d(opt,rcv,n);

	nv=CLK_NUM(opt)*(dck->nobs+1); nc=dck->nx; tnv=nv;
	ix=imat(nx,1);xp=mat(nx,1);Pp=mat(nx,nx);vl=mat(nv,1);vp=mat(nc,1);
	HL=mat(nx,nv);HP=mat(nx,nc);DL=mat(nv,nv);DP=mat(nc,nc);
	for (i=0;i<n;i++) {

		matcpy(xp,dck->x,nx,1);
		matcpy(Pp,dck->P,nx,nx);
		imatcpy(SSC,dck->SSC,MAXSAT,1);
		imatcpy(SRC,dck->SRC,MAXRCV,1);

		/* exclude invalid rcv, unnecessary anymore */
		if (!excInvalidRcv_d(rcv,n,ref)) return 0; 

		/* set or update datum */
		setDatum(dck,rcv,n,flag,SSC,SRC);

		if (!(nv=dckRes(0,dck,rcv,n,nav,xp,Pp,ix,vl,HL,DL,vp,HP,DP,&nc,tnv,flag,NULL,NULL))) {
			trace(2,"%s dck (%d) no valid obs data\n",dck->time,i+1);
			return 0;
		}
		/* lsq filter */
		if (lsqBlock(xp,Pp,ix,vl,HL,DL,vp,HP,DP,nx,nv,nc)) {
			trace(2,"lsq error\n");
			return 0;
		}
		/* get increment */
		for (j=0;j<nx;j++) if (ix[j]) xp[j]+=dck->x[j];

		if (dckRes(i+1,dck,rcv,n,nav,xp,Pp,ix,vl,HL,DL,vp,HP,DP,&nc,tnv,flag,outr,outs)) {
			stat=SOLQ_PPP;
			matcpy(dck->x,xp,nx,1);
			matcpy(dck->P,Pp,nx,nx);
			imatcpy(dck->SSC,SSC,MAXSAT,1);
			imatcpy(dck->SRC,SRC,MAXRCV,1);
			break;
		}
		netRejRcvSat(dck->time,rcv,n,outr,outs);
	}
	if (i>=n) {
		trace(2,"%s dck (%d) iteration overflows\n",dck->time,i);
		states=0;
	}

	if (states&&stat==SOLQ_PPP) {
		if (dckNetAr(dck,rcv,n,xp,Pp)&&
			dckRes(n,dck,rcv,n,nav,xp,Pp,ix,vl,HL,DL,vp,HP,DP,&nc,tnv,flag,outr,
				outs)) {
			ar=1;	/* count fixed ambiguities in current epoch */
		}
		/* track but not fixed solution though ar ok */
		udTrackAmb_d(dck,rcv,n,ar);
		/* fixed solution */
		if (dck->fixedrec>(n/2)) {
			stat=SOLQ_FIX;
			if (ar) {
				matcpy(dck->x,xp,nx,1);
				matcpy(dck->P,Pp,nx,nx);
			}
		}
		else trace(2,"%s not enough fixed rcvs num=%d\n",dck->time,
			dck->fixedrec);

		udSolution_d(dck,rcv,n,stat);
	}

	free(ix);free(xp);free(Pp);
	free(vl);free(HL);free(DL);
	free(vp);free(HP);free(DP);
	return states;
}