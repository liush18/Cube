
#include "dckp.h"

#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

static FILE *fpRclk_dp=NULL;    /* file pointer of ztd */

/* open clk file -------------------------------------------------------------*/
extern int logRclkOpen_dp(const char *file, const prcopt_t *opt)
{
	createdir(file);
	if (!(fpRclk_dp=fopen(file,"w"))) {
		trace(2,"%s: file open error path=%s\n",__func__,file);
		return 0;
	}
	fprintf(fpRclk_dp,"%-22s%6s","% TIME","SOL");
	fprintf(fpRclk_dp,"%20s","RCLK_CODE(s)");
	fprintf(fpRclk_dp,"%20s","RCLK_PHASE(s)");
	fprintf(fpRclk_dp,"%20s","RCLK_WL(c)");
	fprintf(fpRclk_dp,"%20s","CODE_STD");
	fprintf(fpRclk_dp,"%20s","PHASE_STD");
	fprintf(fpRclk_dp,"%20s","WL_STD");

	fprintf(fpRclk_dp,"%20s","RCLK_CODE_FIX(s)");
	fprintf(fpRclk_dp,"%20s","RCLK_PHASE_FIX(s)");
	fprintf(fpRclk_dp,"%20s","RCLK_WL_FIX(c)");
	fprintf(fpRclk_dp,"%20s","CODE_STD_FIX");
	fprintf(fpRclk_dp,"%20s","PHASE_STD_FIX");
	fprintf(fpRclk_dp,"%20s","WL_STD_FIX");
	fprintf(fpRclk_dp,"\n");
	fflush(fpRclk_dp);

	return 1;
}
/* close clk file ------------------------------------------------------------*/
extern void logRclkClose_dp(void)
{
	if (fpRclk_dp) fclose(fpRclk_dp);
	fpRclk_dp=NULL;
}
/* output clk ----------------------------------------------------------------*/
extern void logRclk_dp(rtk_t *rtk)
{
	const prcopt_t *opt=&rtk->opt;
	double *x=rtk->x,*P=rtk->P,*xa=rtk->xa,*Pa=rtk->Pa,fact;
	int i,j,nx=rtk->nx,stat=rtk->sol.stat;

	if (!fpRclk_dp) return;

	fprintf(fpRclk_dp,"%22s%6d",rtk->cTime,stat);
	for (i=0;i<DPNC(opt);i++) {
		fact=i==2?1.0:1E-9;
		j=DPIC(1,i+1,opt);
		fprintf(fpRclk_dp,"%20.12E",x[j]*fact);
	}
	for (i=0;i<DPNC(opt);i++) {
		fact=i==2?1.0:1E-9;
		j=DPIC(1,i+1,opt);
		fprintf(fpRclk_dp,"%20.12E",SQRT(P[j+j*nx])*fact);
	}
	for (i=0;i<DPNC(opt);i++) {
		fact=i==2?1.0:1E-9;
		j=DPIC(1,i+1,opt);
		fprintf(fpRclk_dp,"%20.12E",stat==1?xa[j]*fact:0.0);
	}
	for (i=0;i<DPNC(opt);i++) {
		fact=i==2?1.0:1E-9;
		j=DPIC(1,i+1,opt);
		fprintf(fpRclk_dp,"%20.12E",stat==1?SQRT(Pa[j+j*nx])*fact:0.0);
	}

	fprintf(fpRclk_dp,"\n");
	fflush(fpRclk_dp);
}
