
#include "uck.h"


static FILE *fpRes=NULL;    /* file pointer of residuals */

/* open residuals file -------------------------------------------------------*/
extern int openResLog_u(const char *file)
{
	createdir(file);
	if (!(fpRes=fopen(file, "w"))) {
		trace(2,"res_open: file open error path=%s\n",file);
		return 0;
	}
	return 1;
}
/* close residuals file ------------------------------------------------------*/
extern void closeResLog_u(void)
{
	if (fpRes) fclose(fpRes);
	fpRes=NULL;
}
/* output residuals ----------------------------------------------------------
 * time rcv sat resp resc (resp resc) el -------------------------------------*/
extern void logRes_u(net_t *uck, rcv_t *rcv, const int n)
{
	ssat_t *ssat;
	char id[8];
	int i,j,k,sat;

	if (!fpRes) return;

	for (i=0;i<n&&i<MAXRCV;i++) {

		if (rcv[i].outc) continue;
		for (j=0;j<rcv[i].nobs&&j<MAXOBS;j++) {
			sat=rcv[i].obs[j].sat;

			/* only valid satellite */
			if (rcv[i].datum[sat-1]<OBS_USE) continue;

			satno2id(sat,id);
			ssat=&rcv[i].ssat[sat-1];

			fprintf(fpRes,"%-25s%5s%4s",uck->time,rcv[i].sta.name,id);
			for (k=0;k<UNF(&uck->opt);k++) {
				fprintf(fpRes,"%12.4f%12.4f",ssat->resp[k],ssat->resc[k]);
			}
			fprintf(fpRes,"%12.4f\n",ssat->azel[1]*R2D);
		}
	}

	fflush(fpRes);
}