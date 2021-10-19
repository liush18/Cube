#include "uck.h"

#define SQR(x)          ((x)*(x))
#define SQRT(x)			((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define ROUND(x)        (int)floor((x)+0.5)
#define CONST_AMB       0.001       /* constraint to fixed ambiguity */

#define CONSFIXEDSOL		-1		/* fixed times to fix solution */
#define MAX_AMBSET_STD		1.5		/* max ambiguity std in ambiguity subset */
#define MIN_SAT_SUBSET		4		/* minimum number of ambiguity subset */

/* set ambiguity transfer datum ------------------------------------------------
* notes : call alfter datum initialization by calcuMeasurement_u 
* ---------------------------------------------------------------------------*/
extern void udAmbDatum_u(const prcopt_t *opt, rcv_t *rcv, int n)
{
	ambt_t *amb;
	int i,j,f,sat,fix;

	/* integer ambiguity transfer */
	for (i=0;i<n&&i<MAXRCV;i++) for (j=0;j<rcv[i].nobs&&j<MAXOBS;j++) {
		sat=rcv[i].obs[j].sat;
		if (rcv[i].datum[sat-1]!=OBS_USE) continue;

		amb=&rcv[i].amb[sat-1];
		for (f=0,fix=1;f<UNF(opt);f++) {
			/* count=0 if cycle slip */
			if (amb->count[f]>CONSFIXED) continue;
			else { fix=0; break; }
		}
		/* transfer ambiguity only all frequencies fixed and no slip */
		if (fix) rcv[i].datum[sat-1]=OBS_FIX;
	}
}
/* update tracked ambiguity --------------------------------------------------*/
extern void udTrackAmb_u(net_t *uck, rcv_t *rcv, int n, int ar)
{
	ambt_t *amb;
	int i,j,f,sat,fixed,fixedsat;

	uck->fixedrec=0;
	for (i=0;i<n&&i<MAXRCV;i++) {
		rcv[i].fixed=0;
		if (rcv[i].nobs==0) continue;
		/* integer ambiguity transfer */
		for (j=fixed=fixedsat=0;j<rcv[i].nobs&&j<MAXOBS;j++) {
			sat=rcv[i].obs[j].sat;
			amb=&rcv[i].amb[sat-1];
			/* consider transfered ambiguity as fixed */
			if (rcv[i].datum[sat-1]==OBS_FIX) {
				fixed|=1;
				fixedsat++;
			}

			if (!ar) continue;
			for (f=0;f<UNF(&uck->opt);f++) {
				/* update tracked times */
				if (amb->flag[f]) { /* AR success in current epoch */
					fixed|=1;
					if (amb->value[f]==amb->track[f]) amb->count[f]++;
					else {
						amb->track[f]=amb->value[f];
						amb->count[f]=1;
					}
				}
			}
		}
		/* fixed rcv */
		if (fixedsat>=MIN_SAT_SUBSET) rcv[i].fixed=1;
		else if (uck->sol.stat==SOLQ_FIX) {
			trace(2,"%s float rcv=%s\n",uck->time,rcv[i].sta.name);
		}
		if (fixed) uck->fixedrec+=1;
	}
	if (uck->fixedrec) {
		trace(2,"%s fixed receiver number=%d ar=%d\n",uck->time,uck->fixedrec,ar);
	}
}
/* fixing solution -----------------------------------------------------------
 * return : (0:error, 1:ok) --------------------------------------------------
 * --------------------------------------------------------------------------*/
static int fixAmbSol(const net_t *uck, const rcv_t *rcv, int n, const int *rcvs,
	int num, double *x, double *P)
{
	const prcopt_t *opt=&uck->opt;
	double *v,*H,*R;
	const int nx=uck->nx;
	int i,j,f,sat,nv,info,stat=1,*index,*fixed;

	index=imat(uck->nobs*UNF(opt),1); fixed=imat(uck->nobs*UNF(opt),1);
	/* get fixed ambiguities to fixing solution */
	for (i=nv=0;i<n&&i<MAXRCV;i++) {
		if (rcv[i].nobs==0||!rcvs[i]) continue;
		for (j=0;j<rcv[i].nobs&&j<MAXOBS;j++) {
			sat=rcv[i].obs[j].sat;
			for (f=0;f<UNF(opt);f++) {
				/* estimated and fixed ambiguities, only more than set times */
				if (rcv[i].amb[sat-1].flag[f]&&
					(rcv[i].amb[sat-1].count[f]>CONSFIXEDSOL)) {
					index[nv]=UIAMB(i+1,sat,f+1,opt);
					fixed[nv]=rcv[i].amb[sat-1].value[f];
					nv++;
				}
			}
		}
	}

#if 0
	if (nv!=num||nv==0) {
		if (nv!=num) trace(2,"fixAmbSol: nv!=num nv=%d num=%d\n",nv,num);
		free(index); free(fixed); 
		return 0;
	}
#else
	if (nv==0) {
		free(index); free(fixed); 
		return 0;
	}
#endif

	trace(2,"fixAmbSol: nv = %d\n",nv);

	v=zeros(nv,1); H=zeros(nx,nv); R=zeros(nv,nv);
	/* integer ambiguity constraint */
	for (i=0;i<nv;i++) {
		j=index[i];
		v[i]=(double)fixed[i]-x[j];
		H[j+i*nx]=1.0;
		R[i+i*nv]=SQR(CONST_AMB);
	}

	/* update states with constraints */
	if ((info=filter(x,P,H,v,R,nx,nv))) {
		trace(2,"%s fix solution filter error info=%d\n",uck->time,info);
		stat=0;
	}

	free(index); free(fixed); free(v); free(H); free(R);
	return stat;
}
/* fixing ambiguity by lambda ------------------------------------------------
* return : (-1:AR error, 0: sat subset less than 4, k:fixed ambiguities )
* ----------------------------------------------------------------------------*/
static int fixAmbLambda(const net_t *uck, rcv_t *rcv, int ircv, 
	const int *sats, const double *x, const double *P)
{
	const prcopt_t *opt=&uck->opt;
	double *a,*Qa,*F,test[2];
	const int nx=uck->nx,n=rcv->nobs;
	int i,j,k,index[MAXOBS*NFREQ]={0},stat=0;

	for (i=0;i<MAXSAT;i++) for (j=0;j<UNF(opt);j++) {
		rcv->amb[i].flag[j]=0;
		rcv->amb[i].value[j]=0;
	}

	/* record ambiguity index to be fixed */
	for (i=k=0;i<n&&i<MAXOBS;i++) {
		if (!sats[i]) continue;
		for (j=0;j<UNF(opt);j++) {
			index[k]=UIAMB(ircv,sats[i],j+1,opt);
			k++;
		}
	}
	if (k<MIN_SAT_SUBSET*UNF(opt)) return 0;

	a=mat(k,1); Qa=mat(k,k); F=mat(k,2);
	/* obtain ambiguity values and co-variances */
	for (i=0;i<k;i++) {
		a[i]=x[index[i]];
		for (j=0;j<k;j++) Qa[j+i*k]=P[index[j]+index[i]*nx];
	}
#if 0
	trace(2,"%s rcv=%s\n",uck->time,rcv->sta.name);
	traceimat(2,sats,1,MAXOBS,3);
	tracemat(2,Qa,k,k,8,3);
#endif
	test[1]=opt->thresar[1];
	/* fixing ambiguity by lambda */
	if (!lambdaAR(k,2,a,Qa,F,test)) {
		rcv->nobs=0;		/* exclude rcv if lambda error */
		trace(2,"%s rcv(%s) lambda error\n",uck->time,rcv->sta.name);
		free(a); free(Qa); free(F);
		return 0;
	}
	else {
		/* varidation by ratio-test */
		if (test[0]<opt->thresar[0]||test[1]<opt->thresar[1]) {
			trace(4,"%s rcv(%s) varidation error : n=%d ratio=%.2f rate=%.2f%%\n",
				uck->time,rcv->sta.name,k/UNF(opt),test[0],test[1]*100.0);
		}
		else {
			stat=1;
			rcv->ratio=(float)test[0];
			rcv->rate =(float)(test[1]*100.0);
			trace(4,"%s rcv(%s) varidation ok    : n=%d ratio=%.2f rate=%.2f%%\n",
				uck->time,rcv->sta.name,k/UNF(opt),rcv->ratio,rcv->rate);
			
			for (i=k=0;i<n&&i<MAXOBS;i++) {
				if (!sats[i]) continue;
				for (j=0;j<UNF(opt);j++) {
					/* AR flag and fixed values */
					rcv->amb[sats[i]-1].flag[j]=1;
					rcv->amb[sats[i]-1].value[j]=ROUND(F[k]);
					k++;
				}
			}
		}
	}

	free(a); free(Qa); free(F);
	return stat?k:-1;
}
/* fixing ambiguities --------------------------------------------------------
 * return : (0:fail, >1:success)
 * ---------------------------------------------------------------------------*/
static int fixAmb(const net_t *uck, rcv_t *rcv, int ircv, int *sats, 
	const double *x, const double *P)
{
	const prcopt_t *opt=&uck->opt;
	double std[MAXOBS][NFREQ]={0},stds[MAXOBS];
	int i,f,k;

	/* record ambiguity stds */
	for (i=0;i<rcv->nobs&&i<MAXOBS;i++) {
		/* skip not in subset */
		if (!sats[i]) continue;
		/* record std of each frequency */
		for (f=0;f<UNF(opt);f++) {
			k=UIAMB(ircv,sats[i],f+1,opt);
			std[i][f]=SQRT(P[k+k*uck->nx]);
		}
	}

	/* partial ambiguity fixing */
	while (1) {
		k=fixAmbLambda(uck,rcv,ircv,sats,x,P);
		if (k<0) { /* AR error and reject an element for new AR */
			/* satellite ambiguity std of all frequencies */
			for (i=0;i<rcv->nobs&&i<MAXOBS;i++) {
				stds[i]=0.0;
				for (f=0;f<UNF(opt);f++) stds[i]+=SQR(std[i][f]);
				stds[i]=SQRT(stds[i]);
			}
			/* satellite ambiguity with maximum std */
			i=d_max(stds,rcv->nobs);
			if (i<0) return 0;
			/* reject satellite ambiguity with maximum std from subset */
			sats[i]=0; 
			for (f=0;f<UNF(opt);f++) std[i][f]=0.0;
		}
		else break;
	}

	return k;
}
/* select satellite subset for AR --------------------------------------------*/
static int selSatSubset_u(const net_t *uck, const rcv_t *rcv, int ircv, int *sats,
	const double *P)
{
	const prcopt_t *opt=&uck->opt;
	double elmask;
	int i,j,k,stdflag,sat,count=0;

	for (i=0;i<MAXOBS;i++) sats[i]=0;
	elmask=uck->opt.elmaskar>0.0?uck->opt.elmaskar:uck->opt.elmin;

	/* count satellite for AR */
	for (i=0;i<rcv->nobs&&i<MAXOBS;i++) {
		sat=rcv->obs[i].sat;
		/* skip satellites rejected or with datum and transfer ambiguities*/
		if (rcv->datum[sat-1]!=OBS_USE) continue;
		/* skip low elevation angle satellites */
		if (rcv->ssat[sat-1].azel[1]<elmask) continue;

		/* skip low precision ambiguities */
		for (j=stdflag=0;j<UNF(opt);j++) {
			k=UIAMB(ircv,sat,j+1,opt);
			if (SQRT(P[k+k*uck->nx])>MAX_AMBSET_STD) {stdflag=1; break;}
		}
		if (stdflag) continue;
		/* record satllite subset */
		sats[i]=sat;
		count++;
	}

	if (count<MIN_SAT_SUBSET) return 0;

	return 1;
}
/* ambiguity resolution ------------------------------------------------------
 * return : 1:ok, 0:error 
 * ---------------------------------------------------------------------------*/
extern int uckAmbResol(const net_t *uck, rcv_t *rcv, int n, double *x, double *P)
{
	int i,j,rcvs[MAXRCV]={0},sats[MAXOBS],valid=0;

	/* AR receiver by receiver */
	for (i=0;i<n&&i<MAXRCV;i++) {
		/* skip invalid rcv or fixed rcv */
		if (rcv[i].nobs==0||rcv[i].fixed) continue;
		if (!selSatSubset_u(uck,&rcv[i],i+1,sats,P)) continue;
		/* fixed ambiguity number */
		if ((j=fixAmb(uck,&rcv[i],i+1,sats,x,P))) {
			rcvs[i]=1;
			valid+=j;
		}
	}

	/* fixing solution for all receivers and modify ambiguity state */
	if (valid) {
		if (fixAmbSol(uck,rcv,n,rcvs,valid,x,P)) {
			return valid;
		}
	}

	/* failed AR */
	return 0;
}