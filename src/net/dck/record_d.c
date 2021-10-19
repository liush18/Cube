/* record last epoch states for constraint -----------------------------------*/

#include "dck.h"

#define T_COL		16		/* total columns of output matrix */
#define T_COL_DEC	8		/* columns of output matrix under decimal point */
#define TIME_INTER	30.0	/* time interval (s) */

#define SQR(x)	((x>0)?(x*x):(-(x*x)))
#define SQRT(x)	((x==0)?(0.0):((x>0)?(sqrt(x)):(-sqrt(-x))))

static FILE *fpr1=NULL;
static FILE *fpr2=NULL;
static FILE *fpg1=NULL;
static FILE *fpg2=NULL;
static char fileInd_[MAXSTRPATH];		/* constraint index file */
static char fileVar_[MAXSTRPATH];		/* constraint co-variance file */

/* record output file path ---------------------------------------------------*/
extern void recordOutFile_d(const char *file)
{
	strcpy(fileInd_,file);strcat(fileInd_,".ind");
	strcpy(fileVar_,file);strcat(fileVar_,".var");
	fprintf(stdout,"out: %s\n",fileInd_);
	fprintf(stdout,"out: %s\n",fileVar_);
}
/* open files to record ------------------------------------------------------*/
extern void openRecord_d(int lastEp)
{
	if (lastEp!=LST_EPOCH) return;
	if (!(fpr1=fopen(fileInd_,"w"))||!(fpr2=fopen(fileVar_,"w"))) {
		trace(2,"open record file error %s %s\n",fileInd_,fileVar_);
		return;
	}
}
/* close recorded files ------------------------------------------------------*/
extern void closeRecord_d(int lastEp)
{
	if (lastEp!=LST_EPOCH) return;
	if (fpr1) fclose(fpr1); fpr1=NULL;
	if (fpr2) fclose(fpr2); fpr2=NULL;
}
/* output record information ---------------------------------------------------
* output record information for the next processing
* args   : int		type		I		record information type, see notes
*		   char		*time		I		time string (like 2021/02/02 00:00:00.00)
*		   char		*rcv		I		rcv name
*		   int		sat			I		satellite number
*		   int		frq			I		frequency number
*		   double	value		I		estimated values
*		   int		lastEp		I		last epoch flag
* return : none
* ----------------------------------------------------------------------------*/
extern void record_d(int type, const char *time, const char *rcv, int sat, 
	int frq, int trp, int datum, double value, int lastEp)
{
	char satid[8];

	if (!fpr1||lastEp!=LST_EPOCH) return;

	if (type==TYPE_TIME) {
		/* #TIME,time */
		fprintf(fpr1,"%s %s\n",CST_STR_TIME,time);
	}
	else if (type==TYPE_TRP) {
		/* #TRP,rcv name,value */
		fprintf(fpr1,"%s %s %d %.6f\n",CST_STR_TRP,rcv,trp,value);
	}
	else if (type==TYPE_AMB) {
		/* #AMB,rcv name,satid,frequency number,datum,value */
		satno2id(sat,satid);
		fprintf(fpr1,"%s %s %s %d %d %.6f\n",
			CST_STR_AMB,rcv,satid,frq,datum,value);
	}
	fflush(fpr1);
}
/* output constraint co-variance ----------------------------------------------
* print matrix to stdout
* args  : double *A     I   matrix A (n xp m)
*         int    nc		I   number of rows and columns of DP
*		  int    lastEp I   last epoch (value=2)
* return: none
* ----------------------------------------------------------------------------*/
extern void recordVar_d(const double *DP, int nc, int lastEp)
{
	int i,j;

	if (!fpr2||lastEp!=LST_EPOCH) return;
	for (i=0;i<nc;i++) {
		for (j=0;j<nc;j++) {
			fprintf(fpr2,"%*.*e",T_COL,T_COL_DEC,DP[i+j*nc]);
		}
		fprintf(fpr2,"\n");
	}
	fflush(fpr2);
}
/* open record file to get constraint ----------------------------------------*/
static int openGetRecord(const char *fileInd, const char *fileVar)
{
	if (!(fpg1=fopen(fileInd,"r"))||!(fpg2=fopen(fileVar,"r"))) {
		return 0;
	}
	return 1;
}
/* close recorded files ------------------------------------------------------*/
static void closeGetRecord(void)
{
	if (fpg1) fclose(fpg1); fpg1=NULL;
	if (fpg2) fclose(fpg2); fpg2=NULL;
}
/* initialize csd_t struct ---------------------------------------------------*/
static void initCsd(csd_t *csd)
{
	csd->type=csd->sat=csd->frq=csd->trp=0;
	csd->rcv[0]='\0';
	csd->value=0.0;
}
/* initialize dtm_t struct ---------------------------------------------------*/
static void initDtm(dtm_t *dtm)
{
	dtm->rcv[0]='\0';
	dtm->sat=dtm->frq=dtm->datum=0;
	dtm->value=0.0;
}
/* get record index ----------------------------------------------------------*/
static int getRecordIndex(net_t *uck)
{
	csd_t *csd;
	dtm_t *dtm;
	double value;
	char buff[1024],rcv[8],satid[4];
	int i=0,frq,datum;

	while (fgets(buff,sizeof(buff),fpg1)) {
		if (i==0) {
			if (str_time(buff,6,19,&uck->ctime)) {
				trace(2,"getRecordIndex: record invalid epoch %s\n",buff);
				return 0;
			}
			i++;
			continue;
		}

		if (uck->nc>=uck->ncmax) {
			uck->ncmax+=1024;
			if (!(csd=(csd_t *)realloc(uck->csd,sizeof(csd_t)*(uck->ncmax)))) {
				trace(1,"getRecordIndex: malloc error ncmax=%d\n",uck->ncmax);
				free(uck->csd); uck->csd=NULL; uck->nc=uck->ncmax=0;
				return 0;
			}
			uck->csd=csd;
		}

		if (uck->nd>=uck->ndmax) {
			uck->ndmax+=1024;
			if (!(dtm=(dtm_t *)realloc(uck->dtm,sizeof(dtm_t)*(uck->ndmax)))) {
				trace(1,"getRecordIndex: malloc error ndmax=%d\n",uck->ndmax);
				free(uck->dtm); uck->dtm=NULL; uck->nd=uck->ndmax=0;
				return 0;
			}
			uck->dtm=dtm;
		}

		initCsd(uck->csd+uck->nc); initDtm(uck->dtm+uck->nd);
		if (!strncmp(CST_STR_TRP,buff,4)) {
			if (sscanf(buff+5,"%s %d %lf",uck->csd[uck->nc].rcv,
				&uck->csd[uck->nc].trp,&uck->csd[uck->nc].value)<2) {
				
				trace(2,"decode %s error %s\n",CST_STR_TRP,buff);
				return 0;
			}
			uck->csd[uck->nc].type=TYPE_TRP;
		}
		else if (!strncmp(CST_STR_AMB,buff,4)) {
			if (sscanf(buff+5,"%s %s %d %d %lf",rcv,satid,&frq,&datum,
				&value)<5) {
				trace(2,"decode %s error %s\n",CST_STR_AMB,buff);
				return 0;
			}
			if (datum>OBS_USE) {
				strncpy(uck->dtm[uck->nd].rcv,rcv,4);
				uck->dtm[uck->nd].rcv[4]='\0';
				uck->dtm[uck->nd].sat=satid2no(satid);
				uck->dtm[uck->nd].frq=frq;
				uck->dtm[uck->nd].datum=datum;
				uck->dtm[uck->nd].value=value;
				uck->nd++;
				continue;
			}
			uck->csd[uck->nc].type=TYPE_AMB;
			strncpy(uck->csd[uck->nc].rcv,rcv,4);
			uck->csd[uck->nc].rcv[4]='\0';
			uck->csd[uck->nc].sat=satid2no(satid);
			uck->csd[uck->nc].frq=frq;
			uck->csd[uck->nc].value=value;
		}
		uck->nc++;
	}

	return 1;
}
/* get record co-variance ----------------------------------------------------*/
static int getRecordVar(net_t *uck)
{
	char *p,*q;
	int i,j,n=T_COL,nc=uck->nc;

	if (!(p=(char*)malloc(sizeof(char)*nc*(n+1)))) {
		trace(2,"getRecordVar: malloc error\n");
		return 0;
	}
	uck->DP=mat(nc,nc);

	i=0;
	while (fgets(p,sizeof(char)*nc*(n+1),fpg2)) {
		q=p;
		if (strlen(p)-1<sizeof(char)*nc*n) {
			trace(2,"getRecordVar: columns error d(var)=%d d(index)=%d\n",
				strlen(p)-1,sizeof(char)*nc*n);
			break;
		}
		for (j=0;j<nc;j++,q+=n) {
			uck->DP[i+j*nc]=(double)str2num(q,0,n);
		}
		i++;
	}
	if (i!=nc) {
		trace(2,"getRecordVar: dimension error d(var)=%d d(index)=%d\n",i,nc);
		for (i=0;i<nc;i++) uck->DP[i]=0.0;
		return 0;
	}
#if 0
	FILE *fp=fopen("dd.txt","w");
	matfprint(uck->DP,nc,nc,m,6,fp);
	fclose(fp);
#endif
	free(p);
	return 1;
}
/* get record information ----------------------------------------------------*/
extern int getRecord_d(net_t *uck, const char *fileInd, const char *fileVar)
{
	if (!openGetRecord(fileInd,fileVar)) {
		trace(2,"getRecord: open file error %s %s\n",fileInd,fileVar);
		closeGetRecord();
		return 0;
	}
	if (!getRecordIndex(uck)||!getRecordVar(uck)) {
		closeGetRecord();
		return 0;
	}
	closeGetRecord();
	return 1;
}
/* copy constraint information to states and co-variances --------------------*/
extern int copyConstraint_d(net_t *uck, rcv_t *rcv, int n)
{
	const prcopt_t *opt=&uck->opt;
	csd_t *csd=uck->csd;
	dtm_t *dtm=uck->dtm;
	gtime_t time;
	char sta[MAXRCV*4+1],*p;
	int i,j,r,nx=uck->nx,nc=uck->nc,nd=uck->nd,nm,*ind,crcv[MAXRCV]={0};
	
	time=timeadd(uck->ctime,TIME_INTER);
	if (fabs(timediff(time,uck->sol.time))>DTTOL) return 0;	/* invalid constraint */

	/* receiver name set */
	for (i=0,p=sta;i<n;i++,p+=4) strncpy(p,rcv[i].sta.name,4);
	*p='\0';

	/* mark states with constraint */
	nm=nc>nd?nc:nd;
	ind=imat(nm,1);
	for (i=0;i<nm;i++) ind[i]=-1;

	for (i=0;i<nc;i++) {
		if (csd[i].type==TYPE_TRP) {
			if (!(p=strstr(sta,csd[i].rcv))) continue;
			r=(int)((p-sta)/4+1);
			if (opt->tropopt==TROPOPT_EST&&csd[i].trp>0) continue;
			crcv[r-1]|=1;
			ind[i]=DIT(r,opt)+csd[i].trp;
		}
		else if (csd[i].type==TYPE_AMB) {
			if (!(p=strstr(sta,csd[i].rcv))) continue;
			r=(int)((p-sta)/4+1);
			crcv[r-1]|=1;
			ind[i]=DIAMB(r,csd[i].sat,csd[i].frq,opt);
		}
		else {
			trace(2,"copyConstraint: invalid constraint\n");
			free(ind);
			return 0;
		}
	}
	/* copy constraint to states and co-variances */
	for (i=0;i<nc;i++) {
		if (ind[i]<0) continue;
		uck->x[ind[i]]=csd[i].value;
		for (j=0;j<nc;j++) {
			if (ind[j]<0) continue;
			uck->P[ind[i]+ind[j]*nx]=uck->DP[i+j*nc];
		}
	}
	/* treat receiver without constraint as new rcv */
	for (i=0;i<n&&i<MAXRCV;i++) {
		if (crcv[i]) continue; /* with constraint */
		rcv[i].outc=MAXRCVOUT+1;
	}

	/* set ambiguity datum */
	for (i=0;i<nm;i++) ind[i]=-1;
	for (i=0;i<MAXRCV;i++) crcv[i]=0;
	for (i=0;i<nd;i++) {
		if (!(p=strstr(sta,dtm[i].rcv))) continue;
		r=(int)((p-sta)/4+1);

		/* record ambiguity datum */
		if (dtm[i].datum==OBS_SRC) uck->SRC[r-1]=dtm[i].sat;
		else if (dtm[i].datum==OBS_SSC) uck->SSC[dtm[i].sat-1]=r;
		else if (dtm[i].datum==OBS_FIX) {
			for (j=0;j<DNF(opt);j++) {
				rcv[r-1].amb[dtm[i].sat-1].count[j]=CONSFIXED+1;
			}
			crcv[r-1]|=1;
		}
		else {
			trace(2,"copyConstraint: invalid ambiguity datum\n");
			continue;
		}
		ind[i]=DIAMB(r,dtm[i].sat,dtm[i].frq,opt);
	}

	for (i=0;i<nd;i++) {
		if (ind[i]<0) continue;
		uck->x[ind[i]]=round(dtm[i].value);
	}
	uck->fixedrec=0;
	for (i=0;i<MAXRCV;i++) if (crcv[i]) uck->fixedrec++;
	free(ind);
	return 1;
}
