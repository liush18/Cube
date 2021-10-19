
#include "uck.h"

#define SQR(x)		(x*x)
#define VAR_GRA     SQR(0.01)   /* init variance gradient (m^2) */
#define GAP_RESTROP	120			/* default gap to reset trops parameters (ep) */
#define GAP_RESION  120         /* default gap to reset ionos parameters (ep) */



//
///* temporal update of tropospheric delay -------------------------------------*/
//static updateClk(net_t *uck, rcv_t *rcv, int n, const nav_t *nav)
//{
//	sol_t sol,sol0={{0}};
//	ssat_t ssats[MAXSAT];
//	double dts;
//	char msg[128]="";
//	int i,k;
//
//	for (i=0;i<MAXSAT;i++) {
//		k=UISC(i+1,&uck->opt);
//		if (uck->x[k]!=0.0) continue;
//		if (!ephclk(uck->sol.time,uck->sol.time,i+1,nav,&dts)) {
//			trace(2,"updateClk error\n");
//			continue;
//		}
//		uck->x[k]=dts;
//	}
//	for (i=0;i<n&&i<MAXRCV;i++) {
//		k=(UIRC(i+1,&uck->opt));
//		if (uck->x[k]!=0.0) continue;
//		sol=sol0;
//		if (!pntpos(rcv[i].obs,rcv[i].nobs,nav,&uck->opt,&sol,NULL,ssats,msg)) {
//			trace(2,"%s spp err rcv(%s) msg(%s)\n",uck->time,rcv[i].sta.name,msg);
//			rcv[i].nobs=0;
//			continue;
//		}
//		uck->x[k]=sol.dtr[0];
//	}
//}
/* temporal update of tropospheric delay -------------------------------------*/
static void updateTrop(net_t *uck, const rcv_t *rcv, int n)
{
	const prcopt_t *opt=&uck->opt;
	double pos[3],azel[]={0.0,PI/2.0},ztd,var;
	int i,j,k,nx=uck->nx;

	for (i=0;i<n&&i<MAXRCV;i++) {
		if (rcv[i].nobs==0) continue;
		k=UIT(i+1,&uck->opt);
		if (uck->x[k]==0.0||(int)rcv[i].outc>GAP_RESTROP) {
			ecef2pos(rcv[i].sta.pos,pos);
			ztd=tropMops(rcv[i].time,pos,azel,&var); /* zenith total delay */
			initxNet(uck,ztd,var,k);

			if (opt->tropopt==TROPOPT_ESTG) {
				for (j=k+1;j<k+3;j++) initxNet(uck,1E-6,VAR_GRA,j);
			}
		}
		else {
			uck->P[k*(nx+1)]+=opt->prn[1]*fabs(rcv[i].tt);
			if (opt->tropopt==TROPOPT_ESTG) {
				for (j=k+1;j<k+3;j++) {
					uck->P[j*(nx+1)]+=SQR(0.1)*opt->prn[1]*fabs(rcv[i].tt);
				}
			}
		}
	}
}
/* temporal update of ionospheric delay --------------------------------------*/
static void updateIono(net_t *uck, const rcv_t *rcv, int n)
{
	const prcopt_t *opt=&uck->opt;
	int i,j,k,sat;

	for (i=0;i<n&&i<MAXRCV;i++) for (j=0;j<rcv[i].nobs&&j<MAXOBS;j++) {
		sat=rcv[i].obs[j].sat;
		k=UII(i+1,sat,opt);
		if (uck->x[k]!=0.0&&(int)rcv[i].ssat[sat-1].outc[0]>GAP_RESION) {
			initxNet(uck,0.0,0.0,k);
		}
		else if (uck->x[k]!=0.0) {
			uck->P[k+k*uck->nx]+=opt->prn[2]*fabs(rcv[i].tt);
		}
	}
}
/* update states of ambiguities ----------------------------------------------*/
static void updateAmb(net_t *uck, rcv_t *rcv, int n)
{
	const prcopt_t *opt=&uck->opt;
	int i,j,k,f,sat,slip;

	/* reset ambiguity of outlock receiver in previous epoch */
	for (i=0;i<n&&i<MAXRCV;i++) if (rcv[i].outc>=MAXRCVOUT) {
		for (j=0;j<MAXSAT;j++) for (f=0;f<UNF(opt);f++) {
			initxNet(uck,0.0,0.0,UIAMB(i+1,j+1,f+1,opt));
			rcv[i].amb[j].count[f]=0;
		}
	}

	/* reset ambiguity of outlock satellite in previous epoch */
	for (j=0;j<MAXSAT;j++) if (uck->outc[j]>=MAXSATOUT) {
		for (i=0;i<n&&i<MAXRCV;i++) for (f=0;f<UNF(opt);f++) {
			initxNet(uck,0.0,0.0,UIAMB(i+1,j+1,f+1,opt));
			rcv[i].amb[j].count[f]=0;
		}
	}

	/* reset ambiguity if outlock satellite for receiver in previous epoch */
	for (i=0;i<n&&i<MAXRCV;i++) for (j=0;j<MAXSAT;j++) {
		for (f=0;f<UNF(opt);f++) {
			if (rcv[i].ssat[j].outc[f]>=MAXSATOUT) {
				initxNet(uck,0.0,0.0,UIAMB(i+1,j+1,f+1,opt));
				rcv[i].amb[j].count[f]=0;
			}
		}
	}

	/* reset ambiguity if cycle slip in current epoch */
	for (i=0;i<n&&i<MAXRCV;i++) for (j=0;j<rcv[i].nobs&&j<MAXOBS;j++) {
		sat=rcv[i].obs[j].sat;
		for (f=slip=0;f<UNF(opt);f++) slip|=rcv[i].ssat[sat-1].slip[f];
		
		for (f=0;f<UNF(opt);f++) {
			k=UIAMB(i+1,sat,f+1,opt);
			if (uck->x[k]!=0.0&&!slip) {
				uck->P[k+k*uck->nx]+=opt->prn[3]*fabs(rcv[i].tt);
			}
			else {
				initxNet(uck,0.0,0.0,k);
				rcv[i].amb[sat-1].count[f]=0;
			}
		}
	}
}
/* temporal update of satellite and receiver delay ---------------------------*/
static void updateDelay(net_t *uck, const rcv_t *rcv, int n)
{
	const prcopt_t *opt=&uck->opt;
	int i,f,k;

	/* reset delay if receiver outage */
	for (i=0;i<n&&i<MAXRCV;i++) for (f=0;f<UNF(opt);f++) {
		k=UIRPB(i+1,f+1,opt);
		if (uck->x[k]!=0.0&&(int)rcv[i].outc>=MAXRCVOUT) initxNet(uck,0.0,0.0,k);
		else uck->P[k+k*uck->nx]+=opt->prn[4]*fabs(rcv[i].tt);
	}

	/* reset delay if satellite outage */
	for (i=0;i<MAXSAT;i++) for (f=0;f<UNF(opt);f++) {
		k=UISPB(i+1,f+1,opt);
		if (uck->x[k]!=0.0&&uck->outc[i]>=MAXSATOUT) initxNet(uck,0.0,0.0,k);
		else uck->P[k+k*uck->nx]+=opt->prn[5]*fabs(rcv[0].tt);
	}
}
/* temporal update of satellite and receiver code bias -----------------------*/
static void updateCob(net_t *uck, const rcv_t *rcv, int n)
{
	const prcopt_t *opt=&uck->opt;
	int i,f,k;

	/* reset delay if receiver outage */
	if (opt->rcb)
	for (i=0;i<n&&i<MAXRCV;i++) for (f=0;f<UNF(opt);f++) {
		k=UIRCB(i+1,f+1,opt);
		if (uck->x[k]!=0.0&&(int)rcv[i].outc>=MAXRCVOUT) initxNet(uck,0.0,0.0,k);
		else uck->P[k+k*uck->nx]+=opt->prn[6]*fabs(rcv[i].tt);
	}

	/* reset delay if satellite outage */
	if (opt->scb)
	for (i=0;i<MAXSAT;i++) for (f=0;f<UNF(opt);f++) {
		k=UISCB(i+1,f+1,opt);
		if (uck->x[k]!=0.0&&uck->outc[i]>=MAXSATOUT) initxNet(uck,0.0,0.0,k);
		else uck->P[k+k*uck->nx]+=opt->prn[7]*fabs(rcv[0].tt);
	}
}
/* update decoupled model states ---------------------------------------------*/
extern void udStates_u(net_t *uck, rcv_t *rcv, int n, const nav_t *nav)
{
	trace(3,"udStates_u\n");

	//updateClk(uck,rcv,n,nav);
	updateTrop (uck,rcv,n);
	updateAmb  (uck,rcv,n);
	if (uck->opt.ionoopt==IONOOPT_EST) {
		updateIono (uck,rcv,n);
		updateDelay(uck,rcv,n);
	}
	updateCob(uck,rcv,n);
}