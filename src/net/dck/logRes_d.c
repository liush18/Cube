
#include "dck.h"

static FILE *fpRes_d=NULL;    /* file pointer of residuals */

#define PRINT_RES_VAR	0

/* open residuals file -------------------------------------------------------*/
extern int openResLog_d(const char *file)
{
	createdir(file);
	if (!(fpRes_d=fopen(file,"w"))) {
		trace(2,"res_open: file open error path=%s\n",file);
		return 0;
	}
	fprintf(fpRes_d,"%-22s","% TIME(GPST)");
	fprintf(fpRes_d,"%5s",  "RCV");
	fprintf(fpRes_d,"%5s",  "SAT");
	fprintf(fpRes_d,"%14s", "CODE");
	fprintf(fpRes_d,"%14s", "PHASE");
	fprintf(fpRes_d,"%14s", "WIDE_LANE");
#if PRINT_RES_VAR
	fprintf(fpRes_d,"%9s", "STD_C");
	fprintf(fpRes_d,"%9s", "STD_P");
	fprintf(fpRes_d,"%9s", "STD_WL");
#endif
	fprintf(fpRes_d,"%9s", "AZ");
	fprintf(fpRes_d,"%9s", "EL");
	fprintf(fpRes_d,"%6s", "DTM");
	fprintf(fpRes_d,"%6s", "OUTC");
	fprintf(fpRes_d,"\n");

	return 1;
}
/* close residuals file ------------------------------------------------------*/
extern void closeResLog_d(void)
{
	if (fpRes_d) fclose(fpRes_d);
	fpRes_d=NULL;
}
/* output residuals ----------------------------------------------------------*/
extern int logRes_d(net_t *dck, rcv_t *rcv, const int n)
{
	ssat_t *ssat;
	char id[4];
	int i,j,sat;

	for (i=0;i<n&&i<MAXRCV;i++) {
		if (rcv[i].outc) continue;
		for (j=0;j<rcv[i].nobs&&j<MAXOBS;j++) {

			sat=rcv[i].obs[j].sat;
			satno2id(sat,id);
			ssat=&rcv[i].ssat[sat-1];

			fprintf(fpRes_d,"%22s",dck->time);
			fprintf(fpRes_d,"%5s",rcv[i].sta.name);
			fprintf(fpRes_d,"%5s",id);
			fprintf(fpRes_d,"%14.4f",ssat->resp[0]);
			fprintf(fpRes_d,"%14.4f",ssat->resc[0]);
			fprintf(fpRes_d,"%14.4f",ssat->resw[0]);
#if PRINT_RES_VAR
			fprintf(fpRes_d,"%9.4f",sqrt(ssat->resp[1]));
			fprintf(fpRes_d,"%9.4f",sqrt(ssat->resc[1]));
			fprintf(fpRes_d,"%9.4f",sqrt(ssat->resw[1]));
#endif
			fprintf(fpRes_d,"%9.2f",ssat->azel[0]*R2D);
			fprintf(fpRes_d,"%9.2f",ssat->azel[1]*R2D);
			fprintf(fpRes_d,"%6d",rcv[i].datum[sat-1]);
			fprintf(fpRes_d,"%6d\n",rcv[i].ssat[sat-1].outc[0]);
		}
	}

	fflush(fpRes_d);
	return 0;
}