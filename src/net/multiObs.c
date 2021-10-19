
#include "net.h"

static FILE *fps_obs[MAXRCV]={NULL};

#define MAXLENPATH		1024

static pcvs_t pcv_s={0};		/* receiver antenna parameters */
static obs_t  obs_s={0};		/* observation data */
static nav_t  nav_s={0};		/* navigation data */
static rcv_t  rcv_s[MAXRCV];	/* receiver(station) information */

extern pcvs_t *getpcvs_t(){ return &pcv_s; }
extern obs_t *getobs_t()  { return &obs_s; }
extern nav_t *getnav_t()  { return &nav_s; }
extern rcv_t *getrcv_t()  { return rcv_s; }

/* read observation file path --------------------------------------------------
* args   : char *file   I       observation file path list
*		   FILE **fps	O		multi-receiver observation file pointers
* return : number of file pointers
* notes  : files and fps one to one correspondence, valid or invalid
*-----------------------------------------------------------------------------*/
extern int openMultiOBSFile(const char *file)
{
	FILE *fp;
	char buff[MAXLENPATH];
	int n=0;

	/* fp=NULL if open error */
	if (!(fp=fopen(file,"r"))) {
		trace(2,"open multi-receiver observation files error: %s\n",file);
		return 0;
	}

	/* open observation files */
	while (fgets(buff,MAXLENPATH,fp)) {
		if (buff[strlen(buff)-1]=='\n') buff[strlen(buff)-1]='\0';
		if (!(fps_obs[n]=fopen(buff,"r"))) {
			trace(2,"openMultiOBSFile: open file error: %s\n",file);
		}
		if (++n>=MAXRCV) break;
	}

	/* close observation file list */
	fclose(fp);
	return n;
}
/* observation RINEX header initialization -----------------------------------*/
static void initRNXHO(rnxo_t *rnx)
{
	int i,j,k;

	rnx->ver=2.10;
	rnx->sys=SYS_GPS;
	rnx->tsys=TSYS_GPS;
	for (i=0;i<NUMSYS;i++) for (j=0;j<MAXOBSTYPE;j++) {
		for (k=0;k<4;k++) rnx->tobs[i][j][k]='\0';
	}
	rnx->flag=0;
}
/* get multi-receiver observation file header ----------------------------------
* args   : FILE  **fps	IO		multi-receiver observation file pointers
*		   rcv_t *rcv	O		receiver strcuts
*		   int   n		I		number of file pointers
* return : none
* notes  : fps and rcv one to one correspondence, fps[i] is closed if invalid 
*		   observation header, and rcv[i].excrcv=1
*-----------------------------------------------------------------------------*/
extern void getMultiOBSHeader(rcv_t *rcv, int n)
{
	FILE *fp;
	rnxo_t *rnx;
	char type=' ';
	int i;

	/* n<MAXRCV, but still adding i<MAXRCV for attention */
	for (i=0;i<n&&i<MAXRCV;i++) {

		fp=fps_obs[i];
		rnx=&rcv[i].rnx;

		/* open error observation file */
		if (fp==NULL) {
			rcv[i].excrcv=1;
			continue;
		}

		initRNXHO(rnx);
		/* close file pointer and set NULL if read header error */
		if (!readrnxh(fp,&rnx->ver,&type,&rnx->sys,&rnx->tsys,rnx->tobs,
			NULL,&rcv[i].sta)) {
			fclose(fp); fp=NULL;
			rcv[i].excrcv=1;
			trace(2,"getMultiOBSHeader: read header error ircv=%d\n",i+1);
			continue;
		}
		rcv[i].sta.name[4]='\0';
	}
}
/* get multi-receiver observation file data body -----------------------------*/
extern int getMultiOBSBody(gtime_t time, const char *opt, rcv_t *rcv, int n)
{
	obsd_t *data;
	rnxo_t *rnx;
	uint8_t slips[MAXSAT][NFREQ+NEXOBS]={{0}};
	int i,j,flag,valid=0;

	for (i=0;i<n&&i<MAXRCV;i++) {

		rcv[i].nobs=0;

		/* skip invalid receiver */
		if (rcv[i].excrcv||fps_obs[i]==NULL) continue;

		data=rcv[i].obs;
		rnx=&rcv[i].rnx;

		/* previous epoch has no obs, but read current or next epoch obs */
		if (fabs(timediff(time,rnx->time))<DTTOL) {
			rcv[i].nobs=rnx->n;
			continue;
		}
		else if (timediff(rnx->time,time)>DTTOL) continue;

		while ((rnx->n=readrnxobsb(fps_obs[i],opt,rnx->ver,&rnx->tsys,rnx->tobs,
			&rnx->flag,data,&rcv[i].sta))) {

			/* end of current file pointer */
			if (rnx->n<0) {
				rnx->n=0;
				fclose(fps_obs[i]);
				fps_obs[i]=NULL;
				break;
			}
			
			flag=0;
			for (j=0;j<rnx->n;j++) {
				if (rnx->tsys==TSYS_UTC) data[j].time=utc2gpst(data[j].time);
				saveslips(slips,data+j);
			}

			if (timediff(time,data[0].time)>DTTOL) continue;
			else if (timediff(data[0].time,time)>DTTOL) flag=1;

			for (j=0;j<rnx->n;j++) {
				restslips(slips,data+j);
				data[j].rcv=(uint8_t)i;
			}

			if (!flag) rcv[i].nobs=rnx->n;

			rnx->time=data[0].time;
			valid++;
			break;
		}
	}
	if (valid) return 1;
	return 0;
}
/* get multi-receiver observations ---------------------------------------------
* args   : char *file   I		observation file path list
* return : status (1:ok,0:no data,-1:error)
*-----------------------------------------------------------------------------*/
extern int getMultiOBS(rcv_t *rcv, const char *file)
{
	FILE *fp;
	gtime_t ts;
	obsd_t *obs;
	double eps[6]={2020.0,8.0,30.0, 0.0, 0.0, 0.0};
	int tint=30;
	int i,j,n;

	ts=epoch2time(eps);

	/* open multi-receiver observation file pointers */
	if ((n=openMultiOBSFile(file))==0) {
		trace(1,"no multi-receiver observation files: %s\n",file);
		return -1;
	}

	getMultiOBSHeader(rcv,n);

	fp=fopen("C:\\Users\\Lenovo\\Desktop\\obs.txt","w");
	while (getMultiOBSBody(ts,"",rcv,n)) {
		ts.time+=tint;

		for (i=0;i<n;i++) {
			if (rcv[i].nobs==0) continue;
			obs=rcv[i].obs;
			rcv[i].sta.name[4]='\0';
			for (j=0;j<rcv[i].nobs;j++) {
				fprintf(fp,"%s \t %s \t %d \t %f \t %f \t %f \t %f\n",
					time_str(obs[j].time,2),rcv[i].sta.name,
					obs[j].sat,obs[j].P[0],obs[j].L[0],
					obs[j].P[1],obs[j].L[1]);
			}
		}

		fprintf(fp,"\n\n");
	}
	fclose(fp);
	return 1;
}
/* reject special sat observed by rcvs ---------------------------------------*/
extern void netRejSatAllRcv(rcv_t *rcv, int n, int rej_sat)
{
	int i,j,sat;

	for (i=0;i<n&&i<MAXRCV;i++) {
		for (j=0;j<rcv[i].nobs&&j<MAXOBS;j++) {
			sat=rcv[i].obs[j].sat;
			if (sat=rej_sat) rcv[i].datum[sat-1]=OBS_REJ;
		}
	}
}
/* reject all satellite of single receiver -----------------------------------*/
extern void netRejRcvAllSat(rcv_t *rcv)
{
	int i,sat;

	for (i=0;i<rcv->nobs&&i<MAXOBS;i++) {
		sat=rcv->obs[i].sat;
		rcv->datum[sat-1]=OBS_REJ;
	}
}
/* reject receiver and satellite -----------------------------------------------
*  notes: outr -> rejected times of receiver for post-residuals
*		  outs -> rejected times of satellite for post-residuals
* ----------------------------------------------------------------------------*/
extern void netRejRcvSat(const char *time, rcv_t *rcv, int n, int *outr, 
	int *outs)
{
	int i,j,sat,count;

	/* reject receiver in network */
	for (i=count=0;i<n&&i<MAXRCV;i++) {
		if (outr[i]>=MAX_REJTIMES) {
			outr[i]=0; count++;
			/* rej all satellites of this rcv */
			for (j=0;j<rcv[i].nobs&&j<MAXOBS;j++) {
				sat=rcv[i].obs[j].sat;
				rcv[i].datum[sat-1]=OBS_REJ;
			}
		}
	}
	/* reject satellite in network */
	for (j=0;j<MAXSAT;j++) {
		if (outs[j]>=MAX_REJTIMES) {
			outs[j]=0; count++;
			for (i=0;i<n&&i<MAXRCV;i++) rcv[i].datum[j]=OBS_REJ;
		}
	}
	if (count) {
		trace(2,"%s resoling for rejected receiver or satellite\n",time);
	}
}
/* exclude observation data ----------------------------------------------------
* exclude satellites and delete obs contaning 0.0, only L1, L2
* args  :rcv_t	   *rcv	   IO     receiver parameters
*			prcopt_t   *opt    I      processing options
*			int			n      I	  number of receivers
* return:none
*----------------------------------------------------------------------------*/
extern void excludeObs(rcv_t *rcv, const prcopt_t *opt, const int n)
{
	obsd_t *obs;
	int i,j,k,sat;

	trace(3,"obsexclude: station number=%d\n",n);

	if (n<=0) return;

	/* exclude satellites and delete data containing 0.0 */
	for (i=0;i<n&&i<MAXRCV;i++) {

		for (j=0;j<MAXSAT;j++) rcv[i].datum[j]=0;

		/* receiver without obs */
		if (rcv[i].nobs==0) continue;
		/* excluded receiver */
		if (rcv[i].excrcv) { netRejRcvAllSat(&rcv[i]); continue; }

		obs=rcv[i].obs;
		/* exclude invalid obs */
		for (j=k=0;j<rcv[i].nobs&&j<MAXOBS;j++) {
			sat=obs[j].sat;
			if (!(satsys(sat,NULL)&opt->navsys)||opt->exsats[sat-1]==1||
				obs[j].P[0]==0.0||obs[j].P[1]==0.0||
				obs[j].L[0]==0.0||obs[j].L[1]==0.0) {
				rcv[i].datum[sat-1]=OBS_REJ;
			}
			else k++; /* count of valid obs */
		}
		if (k<4) {
			netRejRcvAllSat(&rcv[i]);
			trace(2,"%s exclude rcv (%s) with less than 4 sats\n",
				time_str(obs[0].time,2),rcv[i].sta.name);
		}
	}
}
/* output observation type ---------------------------------------------------*/
extern void outObsType(FILE *fp, const rcv_t *rcv, int n)
{
	const obsd_t *obs;
	char str[32];
	int i;

	for (i=0;i<n&&i<MAXRCV;i++) {
		obs=rcv[i].obs;
		time2str(obs[0].time,str,0);
		fprintf(fp,"%s %s %s %s %s %s\n",str,rcv[i].sta.name,
			code2obs(obs[0].code[0]),code2obs(obs[0].code[1]),
			code2obs(obs[1].code[0]),code2obs(obs[1].code[1]));
	}
}
/* output observation data ---------------------------------------------------*/
extern void outObsData(FILE *fp, const rcv_t *rcv, int n)
{
	const obsd_t *obs;
	char str[32],id[4];
	int i,j;

	for (i=0;i<n&&i<MAXRCV;i++) {
		obs=rcv[i].obs;
		time2str(obs[0].time,str,0);
		for (j=0;j<rcv[i].nobs&&j<MAXOBS;j++) {
			satno2id(obs[j].sat,id);
			fprintf(fp,"%s %s %s %15.3f %15.3f %15.3f %15.3f %s %s\n",
				str,rcv[i].sta.name,id,obs[j].P[0],obs[j].P[1],
				obs[j].L[0],obs[j].L[1],code2obs(obs[j].code[0]),
				code2obs(obs[j].code[1]));
		}
	}
}