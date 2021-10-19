

#include "qc.h"

#define MIN(x,y)    ((x)<(y)?(x):(y))

/* set gf coefficient (m) ----------------------------------------------------*/
static double GFCoef(const double tint)
{
	if      (tint<=   1.0) return 0.05;
	else if (tint<=  20.0) return (0.1/20.0*tint+0.05);
	else if (tint<=  60.0) return 0.15;
	else if (tint<= 100.0) return 0.25;
	else				   return 0.35;
}
/* set mw coefficient (cycle) ------------------------------------------------*/
static double MWCoef(const double tint)
{
	if      (tint<= 1.0) return 2.5;
	else if (tint<=20.0) return (2.5/20.0*tint+2.5);
	else if (tint<=60.0) return 5.0;
	else				 return 7.5;
}
/* L1/L2 geometry-free phase measurement -------------------------------------*/
static double GFMeas(const obsd_t *obs, const nav_t *nav)
{
	double freq1,freq2;

	freq1=sat2freq(obs->sat,obs->code[0],nav);
	freq2=sat2freq(obs->sat,obs->code[1],nav);
	if (freq1==0.0||freq2==0.0||obs->L[0]==0.0||obs->L[1]==0.0) return 0.0;
	return (obs->L[0]/freq1-obs->L[1]/freq2)*CLIGHT;
}
/* L1/L2 wide-lane phase measurement -----------------------------------------*/
static double MWMeas(const obsd_t *obs, const nav_t *nav)
{
	double freq1,freq2;

	freq1=sat2freq(obs->sat,obs->code[0],nav);
	freq2=sat2freq(obs->sat,obs->code[1],nav);

	if (freq1==0.0||freq2==0.0||obs->L[0]==0.0||obs->L[1]==0.0||
		obs->P[0]==0.0||obs->P[1]==0.0) return 0.0;

	return (obs->L[0]-obs->L[1])-
		(freq1-freq2)/CLIGHT*(freq1*obs->P[0]+freq2*obs->P[1])/(freq1+freq2);
}
/* detect cycle slip by LLI --------------------------------------------------*/
static void detCycleSlip_LLI(const prcopt_t *opt, ssat_t *ssat, const obsd_t *obs, 
	const int n, const char *rcv)
{
    int i,j;
    
    trace(4,"detslp_ll: n=%d\n",n);
    
    for (i=0;i<n&&i<MAXOBS;i++) for (j=0;j<NFREQ;j++) {
        if (obs[i].L[j]==0.0||!(obs[i].LLI[j]&3)) continue;
        ssat[obs[i].sat-1].slip[j]=1;
        trace(4,"%s detCycleSlip_LLI: (%s) detected sat=%2d freq=%d\n",
			time_str(obs[0].time,2),rcv,obs[i].sat,j+1);
    }
}
/* detect cycle slip by geometry free phase jump -----------------------------*/
static void detCycleSlip_GF(const prcopt_t *opt, ssat_t *ssat, const obsd_t *obs, 
	const int n, const nav_t *nav, const char *rcv, const double gfcoef)
{
	double g0,g1,elev,thres;
	int i,j,sat;

	trace(4,"detCycleSlip_GF:n=%d\n",n);

	for (i=0;i<n&&i<MAXOBS;i++){

		if ((g1=GFMeas(obs+i,nav))==0.0) continue;
		sat=obs[i].sat;
		if ((elev=ssat[sat-1].azel[1]*R2D)<opt->elmin*R2D) continue;
		thres=elev>15.0?gfcoef:(-elev/15.0+2)*gfcoef;
		thres=MIN(thres,1.5);
		
		g0=ssat[sat-1].gf;
		if (g0!=0.0&&fabs(g1-g0)>thres)
		{
			for (j=0;j<NFREQ;j++) ssat[sat-1].slip[j]|=1;
			trace(4,"%s detCycleSlip_GF:(%s) detected sat=%2d elev=%5.2f gf0=%8.3f gf=%8.3f thres=%5.3f\n",
				time_str(obs[0].time,2),rcv,obs[i].sat,elev,g0,g1,thres);
		}
		ssat[obs[i].sat-1].gf=g1;
	}
}
/* detect cycle slip by widelane jump ----------------------------------------*/
static void detCycleSlip_MW(const prcopt_t *opt, ssat_t *ssat, const obsd_t *obs, 
	const int n, const nav_t *nav, const char *rcv, const double mwcoef)
{
	double w0,w1,elev,thres,MIN_ARC_GAP=300.0;
	int i,j,sat;

	trace(4, "detCycleSlip_MW:n=%d\n",n);

	for (i=0;i<n&&i<MAXOBS;i++) {

		if ((w1=MWMeas(obs+i,nav))==0.0) continue;
		sat=obs[i].sat;
		if ((elev=ssat[sat-1].azel[1]*R2D)<opt->elmin*R2D) continue;
		thres=elev>20.0?mwcoef:(-0.1*elev+3)*mwcoef;
		thres=MIN(thres,6.0);

		w0=ssat[sat-1].mw;
		if (w0!=0.0&&fabs(w1-w0)>thres) {
			for (j=0;j<NFREQ;j++) ssat[sat-1].slip[j]|=1;
			trace(4,"%s detCycleSlip_MW:(%s) detected sat=%2d elev=%5.2f mw0=%8.3f mw=%8.3f thres=%5.3f\n",
				time_str(obs[0].time,2),rcv,obs[i].sat,elev,w0,w1,thres);
		}

		/* if slip, reset mw meas */
		if(ssat[sat-1].slip[0]||ssat[sat-1].slip[1]||
			fabs(timediff(ssat[sat-1].pt[0][0],obs[i].time))>MIN_ARC_GAP) {
			ssat[sat-1].mw=0.0;
			ssat[sat-1].imw=0;
		}
		if(ssat[sat-1].imw==0) {
			ssat[sat-1].mw=w1;
			ssat[sat-1].imw++;
			ssat[sat-1].pt[0][0]=obs[i].time;
		}
		else {
			j=ssat[sat-1].imw;
			ssat[sat-1].mw=(w0*j+w1)/((double)j+1);
			ssat[sat-1].imw++;
			ssat[sat-1].pt[0][0]=obs[i].time;
		}
	}
}
/* detect cycle slip ---------------------------------------------------------*/
extern void detCycleSlip(const prcopt_t *opt, ssat_t *ssat, obsd_t *obs, 
	const int n, const nav_t *nav, const char *rcv, const double tt, int *jumpc)
{
	int i,j;
	
	for (i=0;i<MAXSAT;i++) for (j=0;j<NFREQ;j++) ssat[i].slip[j]=0;
	
	/* detect cycle slip by LLI */
	//detCycleSlip_LLI(opt,ssat,obs,n,rcv);
	
	/* detect cycle slip by geometry-free phase jump */
	detCycleSlip_GF(opt,ssat,obs,n,nav,rcv,GFCoef(tt));
	
	/* detect clock jump */
	detClockJump(ssat,jumpc,obs,n,nav,rcv);
	
	/* detect slip by Melbourne-Wubbena linear combination jump */
	detCycleSlip_MW(opt,ssat,obs,n,nav,rcv,MWCoef(tt));
}