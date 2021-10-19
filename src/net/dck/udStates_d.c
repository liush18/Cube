
#include "dck.h"

#define SQR(x)		(x*x)
#define VAR_GRA     SQR(0.01)   /* init variance gradient (m^2) */
#define GAP_RESTROP	120			/* default gap to reset trops parameters (ep) */

/* temporal update of tropospheric delay -------------------------------------*/
static void updateTrop_d(net_t *dck, const rcv_t *rcv, int n)
{
	const prcopt_t *opt=&dck->opt;
	double pos[3],azel[]={0.0,PI/2.0},ztd,var;
	int i,j,k,nx=dck->nx;

	for (i=0;i<n&&i<MAXRCV;i++) {
		if (rcv[i].nobs==0) continue;
		k=DIT(i+1,&dck->opt);
		if (dck->x[k]==0.0||(int)rcv[i].outc>GAP_RESTROP) {
			ecef2pos(rcv[i].sta.pos,pos);
			ztd=tropMops(rcv[i].time,pos,azel,&var); /* zenith total delay */
			initxNet(dck,ztd,var,k);

			for (j=k+1;j<k+DNT(opt);j++) initxNet(dck,1E-6,VAR_GRA,j);
		}
		else {
			dck->P[k*(nx+1)]+=opt->prn[1]*fabs(rcv[i].tt);
			for (j=k+1;j<k+DNT(opt);j++) {
				dck->P[j*(nx+1)]+=SQR(0.1)*opt->prn[1]*fabs(rcv[i].tt);
			}	
		}
	}
}
/* update states of ambiguities ------------------------------------------------
 * notes: once receiver outlock in the network, reset
 *		  once satellite out lock in the network, reset
 * ---------------------------------------------------------------------------*/
static void updateAmb_d(net_t *dck, rcv_t *rcv, int n)
{
	const prcopt_t *opt=&dck->opt;
	const int nx=dck->nx;
	int i,j,k,f,sat,slip;

	/* reset ambiguity of outlock receiver in previous epoch */
	for (i=0;i<n&&i<MAXRCV;i++) if (rcv[i].outc) {
		for (j=0;j<MAXSAT;j++) for (f=0;f<DNF(opt);f++) {
			initxNet(dck,0.0,0.0,DIAMB(i+1,j+1,f+1,opt));
			rcv[i].amb[j].count[f]=0;
		}
	}

	/* reset ambiguity of outlock satellite in previous epoch */
	for (j=0;j<MAXSAT;j++) if (dck->outc[j]) {
		for (i=0;i<n&&i<MAXRCV;i++) for (f=0;f<DNF(opt);f++) {
			initxNet(dck,0.0,0.0,DIAMB(i+1,j+1,f+1,opt));
			rcv[i].amb[j].count[f]=0;
		}
	}

	/* reset ambiguity if outlock satellite of receiver in previous epoch */
	for (i=0;i<n&&i<MAXRCV;i++) for (j=0;j<MAXSAT;j++) {
		for (f=0;f<DNF(opt);f++) {
			if (rcv[i].ssat[j].outc[f]>=MAXSATOUT) {
				initxNet(dck,0.0,0.0,DIAMB(i+1,j+1,f+1,opt));
				rcv[i].amb[j].count[f]=0;
			}
		}
	}

	/* reset ambiguity if cycle slip in current epoch */
	for (i=0;i<n&&i<MAXRCV;i++) for (j=0;j<rcv[i].nobs&&j<MAXOBS;j++) {
		sat=rcv[i].obs[j].sat;
		slip=rcv[i].ssat[sat-1].slip[0]||rcv[i].ssat[sat-1].slip[1];

		for (f=0;f<DNF(opt);f++) {
			k=DIAMB(i+1,sat,f+1,opt);
			if (dck->x[k]!=0.0&&!slip) {
				dck->P[k+k*nx]+=opt->prn[3]*fabs(rcv[i].tt);
			}
			else {
				initxNet(dck,0.0,0.0,k);
				rcv[i].amb[sat-1].count[f]=0;
			}
		}
	}
}
/* update decoupled model states ---------------------------------------------*/
extern void udStates_d(net_t *dck, rcv_t *rcv, int n)
{
	trace(3,"udStates_d\n");

	updateTrop_d(dck,rcv,n);
	updateAmb_d (dck,rcv,n);
}