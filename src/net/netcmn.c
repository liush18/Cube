
#include "net.h"


#define UCKVERSION	"3.00"
#define PROGRAM		"CUBE"
#define RUNBY		"APM"
#define	TIMEZONE	"UTC"
#define CLKSYS		"G"
#define ANTEX		"igs14.atx"
#define A_CENTER		"APM  APM @ CAS"
#define CLKREF		1
#define TRF			"IGb14"
#define TIMESYS		"GPS"


/* initialize state and covariance -------------------------------------------*/
extern void initxNet(net_t *net, double xi, double var, int i)
{
	int j;
	net->x[i]=xi;
	for (j=0;j<net->nx;j++) {
		net->P[i+j*net->nx]=net->P[j+i*net->nx]=i==j?var:0.0;
	}
}
/* get coordinate for stations -----------------------------------------------*/
extern int readSnxNet(const char *file, rcv_t *rcv, const int n)
{
	FILE *fp;
	char buff[1024],stalist[MAXRCV*5+1],code[5],*p;
	int i,flag=0;

	trace(3,"readsnx:file=%s\n",file);

	if (!(fp=fopen(file,"r"))) {
		trace(2,"snx file open error: %s\n",file);
		return 0;
	}

	p=stalist;
	for (i=0;i<n&&i<MAXRCV;i++) {
		if (rcv[i].excrcv) continue;
		rcv[i].sta.pos[0]=rcv[i].sta.pos[1]=rcv[i].sta.pos[2]=0.0;
		strup(rcv[i].sta.name);
		strncpy(p,rcv[i].sta.name,4);
		p+=4;
		sprintf(p++,"%c",' ');
	}

	while (fgets(buff,sizeof(buff),fp)) {
		if (strstr(buff,"+SOLUTION/ESTIMATE")) {
			flag=1;
			continue;
		}
		if (strstr(buff,"-SOLUTION/ESTIMATE")) break;

		strncpy(code,&buff[14],4);code[4]='\0';
		if (flag&&strstr(stalist,code)) {
			for (i=0;i<n&&i<MAXRCV;i++) {
				if (rcv[i].excrcv) continue;
				if (!strncmp(code,rcv[i].sta.name,4)&&!strncmp(&buff[7],"STAX",4)) {
					rcv[i].sta.pos[0]=str2num(buff,47,21);
					break;
				}	
				else if (!strncmp(code,rcv[i].sta.name,4)&&!strncmp(&buff[7],"STAY",4)) {
					rcv[i].sta.pos[1]=str2num(buff,47,21);
					break;
				}	
				else if (!strncmp(code,rcv[i].sta.name,4)&&!strncmp(&buff[7],"STAZ",4)) {
					rcv[i].sta.pos[2]=str2num(buff,47,21);
					break;
				}	
			}

		}
	}

	for (i=0;i<n&&i<MAXRCV;i++) {
		if (rcv[i].excrcv) continue;
		if (rcv[i].sta.pos[0]*rcv[i].sta.pos[1]*rcv[i].sta.pos[2]==0.0) {
			rcv[i].excrcv=1;
			trace(2,"no snx coordinates: sta=%s\n",rcv[i].sta.name);
		}	
	}

	fclose(fp);
	return 1;
}
/* put antenna parameters into nav -------------------------------------------*/
extern void setPcvNet(gtime_t time, prcopt_t *popt, nav_t *nav, 
	const pcvs_t *pcvs, rcv_t *rcv, const int n)
{
	pcv_t *pcv;
	double pos[3],del[3];
	char id[64];
	int i,j;

	/* set satellite antenna parameters */
	for (i=0;i<MAXSAT;i++) {
		if (!(satsys(i+1,NULL)&popt->navsys)) continue;
		if (!(pcv=searchpcv(i+1,"",time,pcvs))) {
			satno2id(i+1,id);
			trace(3,"no satellite antenna pcv: %s\n",id);
			continue;
		}
		nav->pcvs[i]=*pcv;
	}
	for (i=0;i<n&&i<MAXRCV;i++) {
		/* auto set station parameters, input parameters are no longer received */
		if (rcv[i].sta.deltype==1) { /* xyz */
			if (norm(rcv[i].sta.pos,3)>0.0) {
				ecef2pos(rcv[i].sta.pos,pos);
				ecef2enu(pos,rcv[i].sta.del,del);
				for (j=0;j<3;j++) rcv[i].sta.del[j]=del[j];
				rcv[i].sta.deltype=0;
			}
		}
		if (!(pcv=searchpcv(0,rcv[i].sta.antdes,time,pcvs))) {
			trace(2, "no receiver antenna pcv: %s\n",rcv[i].sta.antdes);
			rcv[i].sta.antdes[0]='\0';
			rcv[i].excrcv=1;
			continue;
		}
		strcpy(rcv[i].sta.antdes,pcv->type);
		rcv[i].pcvr=*pcv;
	}
}
/* read ocean tide loading parameters ----------------------------------------*/
extern void readOtlNet(const char *file, rcv_t *rcv, const int n)
{
	for (int i=0;i<n&&i<MAXRCV;i++)
		if (!readblq(file,rcv[i].sta.name,rcv[i].odisp)){
			//rcv[i].excrcv=1;
		}	
}
/* output clk file header (like clk file) ------------------------------------*/
static int outClkHeader2(const prcopt_t *opt, const rcv_t *rcv, int n, 
	const char *atxf, const char *outf, const int *sats)
{
	FILE *fp;
	gtime_t t={0};
	double ep[6];
	char label[32]={0},str[32]={0},*p;
	int i,j,satnum=0;

	if(!(fp=fopen(outf,"w"))){
		trace(2,"%s : open outf error file=%s\n",__func__,outf);
		return 0;
	}

	/* RINEX VERSION / TYPE */
	sprintf(label,"RINEX VERSION / TYPE");
	fprintf(fp,"%9s%11s%s%19s%20s%s\n","3.00","","C","","",label);

	/* PGM / RUN BY / DATE */
	sprintf(label,"PGM / RUN BY / DATE");
	t=timeget(); time2epoch(t,ep);
	sprintf(str,"%4.0f%02.0f%02.0f %02.0f%02.0f%02.0f %s",
		ep[0],ep[1],ep[2],ep[3],ep[4],ep[5],TIMEZONE);
	fprintf(fp,"%-20s%-20s%-20s%s\n",PROGRAM,RUNBY,str,label);

	/* SYS / PCVS APPLIED */
	sprintf(label,"SYS / PCVS APPLIED");
	str[0]='\0';
	if (p=strrchr(atxf,FILEPATHSEP)) {
		strncpy(str,p+1,10);str[10]='\0';
		strup(str);
	}
	fprintf(fp,"%-20s%-40s%s\n",CLKSYS,str,label);

	/* ANALYSIS CENTER */
	sprintf(label,"ANALYSIS CENTER");
	fprintf(fp,"%-60s%s\n",A_CENTER,label);

	/* TIME SYSTEM ID */
	sprintf(label,"TIME SYSTEM ID");
	fprintf(fp,"%6s%54s%s\n",TIMESYS,"",label);

	/* # / TYPES OF DATA */
	sprintf(label,"# / TYPES OF DATA");
	fprintf(fp,"%6d%6s%6s%42s%s\n",2,"AR","AS","",label);

	/* if not zmc but reference rcv datum */
	if (opt->datumType) {
		/* # OF CLK REF */
		sprintf(label,"# OF CLK REF");
		fprintf(fp,"%9d%51s%s\n",CLKREF,"",label);

		/* ANALYSIS CLK REF */
		sprintf(label,"ANALYSIS CLK REF");
		fprintf(fp,"%-5s%-55s%s\n",rcv[0].sta.name,rcv[0].sta.marker,label);
	}

	/* # OF SOLN STA / TRF */
	sprintf(label,"# OF SOLN STA / TRF");
	for (i=j=0;i<n&&i<MAXRCV;i++) if (!rcv[i].excrcv) j++;
	fprintf(fp,"%6d%4s%-50s%s\n",j,"",TRF,label);

	/* SOLN STA NAME / NUM */
	sprintf(label,"SOLN STA NAME / NUM");
	for (i=0;i<n&&i<MAXRCV;i++) {
		if (rcv[i].excrcv) continue;
		fprintf(fp,"%-5s%-20s%11.0f %11.0f %11.0f%s\n",
			rcv[i].sta.name,rcv[i].sta.marker,rcv[i].sta.pos[0]*1E3,
			rcv[i].sta.pos[1]*1E3,rcv[i].sta.pos[2]*1E3,label);
	}

	/* # OF SOLN SATS */
	sprintf(label,"# OF SOLN SATS");
	for (i=0;i<MAXSAT;i++) if (sats[i]) satnum++;
	fprintf(fp,"%6d%54s%s\n",satnum,"",label);

	/* PRN LIST */
	sprintf(label,"PRN LIST");
	for (i=j=0;i<MAXSAT;i++) {
		if (!sats[i]) continue;
		satno2id(i+1,str);
		fprintf(fp,"%-4s",str);
		j++;
		if (j==15) {fprintf(fp,"%s\n",label);j=0;}
	}
	fprintf(fp,"%*s%s\n",(15-j)*4,"",label);

	/* END OF HEADER */
	sprintf(label,"END OF HEADER");
	fprintf(fp,"%60s%s\n","",label);

	fclose(fp);
	return 1;
}
/* get satllites with valid clocks from temporary clock file -----------------*/
static int getClkSats(const char *tclkf, int *sats)
{
	FILE *fp;
	char buff[1024];
	int i,sat,at_least=20;

	if(!(fp=fopen(tclkf,"r"))){
		trace(2,"%s : open temporary clock file error=%s\n",__func__,tclkf);
		return 0;
	}

	for (i=0;i<MAXSAT;i++) sats[i]=0;
	while (fgets(buff,sizeof(buff),fp)) {
		/* skip lines without satellite clock */
		if (strncmp(buff,"AS ",3)) continue;
		if ((sat=satid2no(buff+3))>0) {
			sats[sat-1]++;
		}
	}
	fclose(fp);

	for (i=0;i<MAXSAT;i++) {
		if (sats[i]<at_least) sats[i]=0;
		else sats[i]=1;
	}

	return 1;
}
/* read temporary clock file and generate final clock file -------------------*/
extern int generateFinalClk(const prcopt_t *opt, const rcv_t *rcv, int n, 
	const char *tclkf, const char *atxf, const char *clkf)
{
	FILE *fpr,*fpw;
	char buff[1024];
	int sats[MAXSAT];

	if (!getClkSats(tclkf,sats)) return 0;
	if (!outClkHeader2(opt,rcv,n,atxf,clkf,sats)) return 0;

	if(!(fpr=fopen(tclkf,"r"))){
		trace(2,"%s : read temporary clock file error=%s\n",__func__,tclkf);
		return 0;
	}

	if(!(fpw=fopen(clkf,"a"))){
		trace(2,"%s : open clock file error=%s\n",__func__,clkf);
		return 0;
	}

	while (fgets(buff,sizeof(buff),fpr)) fprintf(fpw,"%s",buff);

	fclose(fpr);fclose(fpw);
	remove(tclkf);

	return 1;
}