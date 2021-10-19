/*------------------------------------------------------------------------------
* detClockJump.c:detect and repair clock jump
*
*          Copyright (C) 2020-2020 by Shuai Liu, All rights reserved.
*
* references :
*    [1] Guo, F., Zhang, X. Real-time clock jump compensation for precise point
*        positioning. GPS Solut 18, 41¨C50 (2014).
*        https://doi.org/10.1007/s10291-012-0307-3
*    [2] Zhou, F., Dong, D., Li, W. et al. GAMP: An open-source software of 
*        multi-GNSS precise point positioning using undifferenced and uncombined
*        observations. GPS Solut 22, 33 (2018). 
*        https://doi.org/10.1007/s10291-018-0699-9
*
* version:$Revision:$ $Date:$
* history:2020/01/28 1.0  new
*-----------------------------------------------------------------------------*/
#include "qc.h"

#define ROUND(x)    (double)(int)floor((x)+0.5)

/* detect and repair clock jump ----------------------------------------------*/
extern void detClockJump(ssat_t *ssat, int *jumpc, obsd_t* obs, const int n, 
	const nav_t* nav, const char *rcv)
{
	double delp,dell,s,k1,k2,ms,freq1,freq2;
	int i,j,sat,vf=0,jf=0,vsat[MAXPRNGPS]={0},nf=2;

	k1=0.001*CLIGHT-3.0*3.0; k2=2.5E-2; s=0.0;

	for (i=0;i<n;i++) {
		sat=obs[i].sat;
		if ((freq1=sat2freq(sat,obs[i].code[0],nav))==0.0) continue;

		/* only detect gps clock jump */
		if (sat>MAXPRNGPS) continue;
		if (obs[i].P[0]*obs[i].P[1]*obs[i].L[0]*obs[i].L[1]==0.0) 
			continue;
		if (ssat[sat-1].ph[0][0]*ssat[sat-1].ph[0][1]*
			ssat[sat-1].pr[0][0]*ssat[sat-1].pr[0][1]==0.0)
			continue;
		if (ssat[sat-1].slip[0]||ssat[sat-1].slip[1]) continue;

		vf++; /* valid flag */
		delp=obs[i].P[0]-ssat[sat-1].pr[0][0];
		dell=(obs[i].L[0]-ssat[sat-1].ph[0][0])*CLIGHT/freq1;

		if (fabs(delp-dell)>k1) { s+=delp-dell; jf++; } /* jump flag */
	}

	if (jf!=0&&jf==vf) {
		ms=1000.0*s/(jf*CLIGHT);
		if (fabs(ms-ROUND(ms))<k2) {
			*jumpc+=(int)ROUND(ms);
			trace(2,"%s detClockJump: (%s) detected jump=%dms\n",
				time_str(obs[0].time,2),rcv,*jumpc);
		}
	}

	for (i=0;i<n;i++) {
		sat=obs[i].sat;
		if ((freq1=sat2freq(sat,obs[i].code[0],nav))==0.0) continue;
		if ((freq2=sat2freq(sat,obs[i].code[1],nav))==0.0) continue;

		if (sat>MAXPRNGPS) continue;
		vsat[sat-1]=1;
		for (j=0;j<nf;j++) {
			ssat[sat-1].ph[0][j]=obs[i].L[j];
			ssat[sat-1].pr[0][j]=obs[i].P[j];
		}

		/* repair clock jump */
		if (obs[i].L[0]!=0.0) obs[i].L[0]+=(*jumpc*freq1)/(1000.0);
		if (obs[i].L[1]!=0.0) obs[i].L[1]+=(*jumpc*freq2)/(1000.0);
	}
	for (i=0;i<MAXPRNGPS;i++) {
		if (vsat[i]!=0) continue; 
		for (j=0;j<nf;j++) {
			ssat[i].ph[0][j]=0.0;
			ssat[i].ph[1][j]=0.0;
			ssat[i].pr[0][j]=0.0;
			ssat[i].pr[1][j]=0.0;
		}
	}
}