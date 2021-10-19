

#include "net/dck/dck.h"
#include "net/uck/uck.h"
#include "user/ppp/ppp.h"
#include "user/dckp/dckp.h"

static char conf[1024]="F:\\conf_net.conf";

/* begin to uck net solving --------------------------------------------------*/
static void processing(int argc,char **argv)
{
	prcopt_t popt;
	solopt_t sopt;
	filopt_t fopt;
	double sec;
	char stime[64],*p;
	int i,week;

	/* load options from configuration file */
	for (i=1;i<argc;i++) {
		if (!strcmp(argv[i],"-conf")&&i+1<argc) {
			sprintf(conf,"%s",argv[++i]);
		}
	}
	/* load configuration options */
	printf("read: %s\n",conf);
	resetsysopts();
	if (!loadopts(conf,sysopts)) return;
	getsysopts(&popt,&sopt,&fopt);

	sec=time2gpst(popt.time,&week);
	for (i=1;i<argc;i++) {
		/* for network */
		if      (!strcmp(argv[i],"-week"))week=atoi(argv[++i]);
		else if (!strcmp(argv[i],"-dow")) sec =(double)atoi(argv[++i])*86400.0;
		else if (!strcmp(argv[i],"-inh")) popt.inh=1;
		else if (!strcmp(argv[i],"-rcb")) popt.rcb=1;
		else if (!strcmp(argv[i],"-scb")) popt.scb=1;
		
		/* for ppp */
		else if (!strcmp(argv[i],"-eratio1")) {
			popt.eratio[0]=atof(argv[++i]);
			printf("\n%f\n",popt.eratio[0]);
		}
		else if (!strcmp(argv[i],"-ts")) {
			p=stime;
			sprintf(p,"%s",argv[++i]);
			p+=strlen(stime);
			sprintf(p," %s",argv[++i]);
			if (str_time(stime,0,64,&popt.ts)<0) return;
		}
		else if (!strcmp(argv[i],"-te")) {
			p=stime;
			sprintf(p,"%s",argv[++i]);
			p+=strlen(stime);
			sprintf(p," %s",argv[++i]);
			if (str_time(stime,0,64,&popt.te)<0) return;
		}
		else if (!strcmp(argv[i],"-obs")) sprintf(fopt.obs,"%s",argv[++i]);

		/* common */
		else if (!strcmp(argv[i],"-clk")) sprintf(fopt.clk,"%s",argv[++i]);
	}
	popt.time=gpst2time(week,sec);

	if (popt.mode==PMODE_NET) {
		printf("processing week(%04d) dow(%d)\n",week,(int)(sec/86400.0));
		if (popt.dck)	dckPost(&popt,&sopt,&fopt);
		else			uckPost(&popt,&sopt,&fopt);
	}
	else if (popt.mode==PMODE_SINGLE) {
		userPostPos_p(popt.ts,popt.te,popt.ti,popt.tu,&popt,&sopt,&fopt);
	}
	else if (popt.mode>=PMODE_PPP_KINEMA||popt.mode<=PMODE_PPP_FIXED) {
		if (popt.dck) {
			userPostPos_dp(popt.ts,popt.te,popt.ti,popt.tu,&popt,&sopt,&fopt);
		}
		else {
			userPostPos_p(popt.ts,popt.te,popt.ti,popt.tu,&popt,&sopt,&fopt);
		}
	}
	else {
		printf("error mode=%d\n",popt.mode);
	}
}
/* main ----------------------------------------------------------------------*/
int main(int argc, char **argv)
{
	long t1,t2;

	t1=clock();
	processing(argc,argv);
	t2=clock();

	printf("\ntotal time: %6.3f seconds\n",(double)(t2-t1)/CLOCKS_PER_SEC);
	if (argc==1) system("pause");
	return 0;
}
/* test main -----------------------------------------------------------------*/
void main2()
{
	setsysopts(NULL,NULL,NULL);
	saveopts("F:\\confprn.conf","w","configure file",sysopts);
	system("pause");
}
void main3()
{
	nav_t nav={0};
	sta_t sta[10]={{0}};
	obs_t obs={0};
	char file[]="C:\\Users\\Lenovo\\Desktop\\uc_test\\*.20o";

	readrnx(file,1,"",&obs,&nav,sta);
}
void main4()
{
	char file[]="C:\\Users\\Lenovo\\Desktop\\uc_test\\*.sp3";
	nav_t nav={0};
	readsp3(file,&nav,0);
}
void main5()
{
	rcv_t rcv[3]={0};
	char file[]="C:\\Users\\Lenovo\\Desktop\\path.txt";
	getMultiOBS(rcv,file);
}
void main6()
{
	printf("%*.*e\n",15,8,0.00000078);
	system("pause");
}
