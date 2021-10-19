
#include "ppp.h"

#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

static FILE *fpRclk_p=NULL;    /* file pointer of ztd */

/* open clk file -------------------------------------------------------------*/
extern int logRclkOpen_p(const char *file, const prcopt_t *opt)
{
	createdir(file);
	if (!(fpRclk_p=fopen(file,"w"))) {
		trace(2,"%s: file open error path=%s\n",__func__,file);
		return 0;
	}
	fprintf(fpRclk_p,"%-22s%6s","% TIME","SOL");
	fprintf(fpRclk_p,"%20s","RCLK(s)");
	fprintf(fpRclk_p,"%20s","RCLK_STD");

	fprintf(fpRclk_p,"%20s","RCLK_FIX(s)");
	fprintf(fpRclk_p,"%20s","RCLK_STD");
	fprintf(fpRclk_p,"\n");
	fflush(fpRclk_p);

	return 1;
}
/* close ztd file ------------------------------------------------------------*/
extern void logRclkClose_p(void)
{
	if (fpRclk_p) fclose(fpRclk_p);
	fpRclk_p=NULL;
}
/* output ztd ----------------------------------------------------------------*/
extern void logRclk_p(rtk_t *rtk)
{
	const prcopt_t *opt=&rtk->opt;
	double *x=rtk->x,*P=rtk->P,*xa=rtk->xa,*Pa=rtk->Pa,fact=1E-9;
	int i,nx=rtk->nx,stat=rtk->sol.stat;

	if (!fpRclk_p) return;

	i=PIC(1,opt);	/* only GPS */
	fprintf(fpRclk_p,"%22s%6d",rtk->cTime,stat);
	fprintf(fpRclk_p,"%20.12E",x[i]*fact);
	fprintf(fpRclk_p,"%20.12E",SQRT(P[i+i*nx])*fact);
	fprintf(fpRclk_p,"%20.12E",stat==1?xa[i]*fact:0.0);
	fprintf(fpRclk_p,"%20.12E",stat==1?SQRT(Pa[i+i*nx])*fact:0.0);
	fprintf(fpRclk_p,"\n");
	fflush(fpRclk_p);
}
