
#include "base.h"


#define SQR(x)      ((x)*(x))
#define MAX_VAR_EPH SQR(300.0)  /* max variance eph to reject satellite (m^2) */


/* satellite system+prn/slot number to satellite number ------------------------
* convert satellite system+prn/slot number to satellite number
* args  :int    sys       I   satellite system (SYS_GPS,SYS_GLO,...)
*          int    prn       I   satellite prn/slot number
* return:satellite number (0:error)
*-----------------------------------------------------------------------------*/
extern int satno(int sys, int prn)
{
	if (prn<=0) return 0;
	switch (sys) {
	case SYS_GPS:
		if (prn<MINPRNGPS||MAXPRNGPS<prn) return 0;
		return prn-MINPRNGPS+1;
	case SYS_GLO:
		if (prn<MINPRNGLO||MAXPRNGLO<prn) return 0;
		return NSATGPS+prn-MINPRNGLO+1;
	case SYS_GAL:
		if (prn<MINPRNGAL||MAXPRNGAL<prn) return 0;
		return NSATGPS+NSATGLO+prn-MINPRNGAL+1;
	case SYS_QZS:
		if (prn<MINPRNQZS||MAXPRNQZS<prn) return 0;
		return NSATGPS+NSATGLO+NSATGAL+prn-MINPRNQZS+1;
	case SYS_CMP:
		if (prn<MINPRNCMP||MAXPRNCMP<prn) return 0;
		return NSATGPS+NSATGLO+NSATGAL+NSATQZS+prn-MINPRNCMP+1;
	case SYS_IRN:
		if (prn<MINPRNIRN||MAXPRNIRN<prn) return 0;
		return NSATGPS+NSATGLO+NSATGAL+NSATQZS+NSATCMP+prn-MINPRNIRN+1;
	case SYS_LEO:
		if (prn<MINPRNLEO||MAXPRNLEO<prn) return 0;
		return NSATGPS+NSATGLO+NSATGAL+NSATQZS+NSATCMP+NSATIRN +
			prn-MINPRNLEO+1;
	case SYS_SBS:
		if (prn<MINPRNSBS||MAXPRNSBS<prn) return 0;
		return NSATGPS+NSATGLO+NSATGAL+NSATQZS+NSATCMP+NSATIRN+NSATLEO +
			prn-MINPRNSBS+1;
	}
	return 0;
}
/* satellite number to satellite system ----------------------------------------
* convert satellite number to satellite system
* args  :int    sat       I   satellite number (1-MAXSAT)
*          int    *prn      IO  satellite prn/slot number (NULL: no output)
* return:satellite system (SYS_GPS,SYS_GLO,...)
*-----------------------------------------------------------------------------*/
extern int satsys(int sat, int *prn)
{
	int sys=SYS_NONE;
	if (sat<=0||MAXSAT<sat) sat=0;
	else if (sat<=NSATGPS) {
		sys=SYS_GPS;sat+=MINPRNGPS-1;
	}
	else if ((sat-=NSATGPS)<=NSATGLO) {
		sys=SYS_GLO;sat+=MINPRNGLO-1;
	}
	else if ((sat-=NSATGLO)<=NSATGAL) {
		sys=SYS_GAL;sat+=MINPRNGAL-1;
	}
	else if ((sat-=NSATGAL)<=NSATQZS) {
		sys=SYS_QZS;sat+=MINPRNQZS-1;
	}
	else if ((sat-=NSATQZS)<=NSATCMP) {
		sys=SYS_CMP;sat+=MINPRNCMP-1;
	}
	else if ((sat-=NSATCMP)<=NSATIRN) {
		sys=SYS_IRN;sat+=MINPRNIRN-1;
	}
	else if ((sat-=NSATIRN)<=NSATLEO) {
		sys=SYS_LEO;sat+=MINPRNLEO-1;
	}
	else if ((sat-=NSATLEO)<=NSATSBS) {
		sys=SYS_SBS;sat+=MINPRNSBS-1;
	}
	else sat=0;
	if (prn) *prn=sat;
	return sys;
}
/* satellite id to satellite number --------------------------------------------
* convert satellite id to satellite number
* args  :char   *id       I   satellite id (nn,Gnn,Rnn,Enn,Jnn,Cnn,Inn or Snn)
* return:satellite number (0: error)
* notes :120-142 and 193-199 are also recognized as sbas and qzss
*-----------------------------------------------------------------------------*/
extern int satid2no(const char *id)
{
	int sys,prn;
	char code;

	if (sscanf(id,"%d",&prn)==1) {
		if (MINPRNGPS<=prn&&prn<=MAXPRNGPS) sys=SYS_GPS;
		else if (MINPRNSBS<=prn&&prn<=MAXPRNSBS) sys=SYS_SBS;
		else if (MINPRNQZS<=prn&&prn<=MAXPRNQZS) sys=SYS_QZS;
		else return 0;
		return satno(sys,prn);
	}
	if (sscanf(id,"%c%d",&code,&prn)<2) return 0;

	switch (code) {
	case 'G': sys=SYS_GPS;prn+=MINPRNGPS-1;break;
	case 'R': sys=SYS_GLO;prn+=MINPRNGLO-1;break;
	case 'E': sys=SYS_GAL;prn+=MINPRNGAL-1;break;
	case 'J': sys=SYS_QZS;prn+=MINPRNQZS-1;break;
	case 'C': sys=SYS_CMP;prn+=MINPRNCMP-1;break;
	case 'I': sys=SYS_IRN;prn+=MINPRNIRN-1;break;
	case 'L': sys=SYS_LEO;prn+=MINPRNLEO-1;break;
	case 'S': sys=SYS_SBS;prn+=100;break;
	default: return 0;
	}
	return satno(sys,prn);
}
/* satellite number to satellite id --------------------------------------------
* convert satellite number to satellite id
* args  :int    sat       I   satellite number
*          char   *id       O   satellite id (Gnn,Rnn,Enn,Jnn,Cnn,Inn or nnn)
* return:none
*-----------------------------------------------------------------------------*/
extern void satno2id(int sat, char *id)
{
	int prn;
	switch (satsys(sat,&prn)) {
	case SYS_GPS: sprintf(id,"G%02d",prn-MINPRNGPS+1);return;
	case SYS_GLO: sprintf(id,"R%02d",prn-MINPRNGLO+1);return;
	case SYS_GAL: sprintf(id,"E%02d",prn-MINPRNGAL+1);return;
	case SYS_QZS: sprintf(id,"J%02d",prn-MINPRNQZS+1);return;
	case SYS_CMP: sprintf(id,"C%02d",prn-MINPRNCMP+1);return;
	case SYS_IRN: sprintf(id,"I%02d",prn-MINPRNIRN+1);return;
	case SYS_LEO: sprintf(id,"L%02d",prn-MINPRNLEO+1);return;
	case SYS_SBS: sprintf(id,"%03d",prn);return;
	}
	strcpy(id,"");
}
/* test excluded satellite -----------------------------------------------------
* test excluded satellite
* args  :int    sat       I   satellite number
*          double var       I   variance of ephemeris (m^2)
*          int    svh       I   sv health flag
*          prcopt_t *opt    I   processing options (NULL: not used)
* return:status (1:excluded,0:not excluded)
*-----------------------------------------------------------------------------*/
extern int satexclude(int sat, double var, int svh, const prcopt_t *opt)
{
	int sys=satsys(sat,NULL);

	if (svh<0) return 1; /* ephemeris unavailable */

	if (opt) {
		if (opt->exsats[sat-1]==1) return 1;	/* excluded satellite */
		if (opt->exsats[sat-1]==2) return 0;	/* included satellite */
		if (!(sys&opt->navsys)) return 1;		/* unselected sat sys */
	}
	if (sys==SYS_QZS) svh&=0xFE;			/* mask QZSS LEX health */
	if (svh) {
		trace(3,"unhealthy satellite: sat=%3d svh=%02X\n",sat,svh);
		return 1;
	}
	if (var>MAX_VAR_EPH) {
		trace(3,"invalid ura satellite: sat=%3d ura=%.2f\n",sat,sqrt(var));
		return 1;
	}
	return 0;
}
/* geometric distance ----------------------------------------------------------
* compute geometric distance and receiver-to-satellite unit vector
* args  :double *rs       I   satellilte position (ecef at transmission) (m)
*          double *rr       I   receiver position (ecef at reception) (m)
*          double *e        O   line-of-sight vector (ecef)
* return:geometric distance (m) (0>:error/no satellite position)
* notes :distance includes sagnac effect (earth rotation effect) correction
*-----------------------------------------------------------------------------*/
extern double geodist(const double *rs, const double *rr, double *e)
{
	double r;
	int i;

	if (norm(rs,3)<RE_WGS84) return -1.0;
	for (i=0;i<3;i++) e[i]=rs[i]-rr[i];
	r=norm(e,3);
	for (i=0;i<3;i++) e[i]/=r;
	return r+OMGE*(rs[0]*rr[1]-rs[1]*rr[0])/CLIGHT;
}
/* relativity correction -----------------------------------------------------*/
extern void relativity(const double *rs, const double *rr, double *r)
{
	*r+=2.0*dot(rs,rs+3,3)/CLIGHT;
}
/* geometric distance variance -------------------------------------------------
* compute geometric distance variance
* args  :double *vars       I   satellilte position variance
*          double *stdr       I   receiver position std (x,y,z) (m)
* return:geometric distance variance
* notes :not an exact value
*-----------------------------------------------------------------------------*/
extern double geodist_var(const double vars, const double *stdr)
{
	return SQR(sqrt(vars)-norm(stdr,3));
}
/* satellite azimuth/elevation angle -------------------------------------------
* compute satellite azimuth/elevation angle
* args  :double *pos      I   geodetic position {lat,lon,h} (rad,m)
*          double *e        I   receiver-to-satellilte unit vevtor (ecef)
*          double *azel     IO  azimuth/elevation {az,el} (rad) (NULL: no output)
*                               (0.0<=azel[0]<2*pi,-pi/2<=azel[1]<=pi/2)
* return:elevation angle (rad)
*-----------------------------------------------------------------------------*/
extern double satazel(const double *pos, const double *e, double *azel)
{
	double az=0.0,el=PI/2.0,enu[3];

	if (pos[2]>-RE_WGS84) {
		ecef2enu(pos,e,enu);
		az=dot(enu,enu,2)<1E-12?0.0:atan2(enu[0],enu[1]);
		if (az<0.0) az+=2*PI;
		el=asin(enu[2]);
	}
	if (azel) {azel[0]=az;azel[1]=el;}
	return el;
}
/* compute dops ----------------------------------------------------------------
* compute DOP (dilution of precision)
* args  :int    ns        I   number of satellites
*          double *azel     I   satellite azimuth/elevation angle (rad)
*          double elmin     I   elevation cutoff angle (rad)
*          double *dop      O   DOPs {GDOP,PDOP,HDOP,VDOP}
* return:none
* notes :dop[0]-[3] return 0 in case of dop computation error
*-----------------------------------------------------------------------------*/
#define SQRT(x)     ((x)<0.0||(x)!=(x)?0.0:sqrt(x))

extern void dops(int ns,const double *azel, double elmin, double *dop)
{
	double H[4*MAXSAT],Q[16],cosel,sinel;
	int i,n;

	for (i=0;i<4;i++) dop[i]=0.0;
	for (i=n=0;i<ns&&i<MAXSAT;i++) {
		if (azel[1+i*2]<elmin||azel[1+i*2]<=0.0) continue;
		cosel=cos(azel[1+i*2]);
		sinel=sin(azel[1+i*2]);
		H[4*n]=cosel*sin(azel[i*2]);
		H[1+4*n]=cosel*cos(azel[i*2]);
		H[2+4*n]=sinel;
		H[3+4*n++]=1.0;
	}
	if (n<4) return;

	matmul("NT",4,4,n,1.0,H,H,0.0,Q);
	if (!matinv(Q,4)) {
		dop[0]=SQRT(Q[0]+Q[5]+Q[10]+Q[15]);	/* GDOP */
		dop[1]=SQRT(Q[0]+Q[5]+Q[10]);		/* PDOP */
		dop[2]=SQRT(Q[0]+Q[5]);				/* HDOP */
		dop[3]=SQRT(Q[10]);					/* VDOP */
	}
}
/* exclude meas of eclipsing satellite (block IIA) ---------------------------*/
extern void testeclipse(const obsd_t *obs, int n, const nav_t *nav, double *rs)
{
	double rsun[3],esun[3],r,ang,erpv[5]={0},cosa;
	int i,j;
	const char *type;

	trace(3,"testeclipse:\n");

	/* unit vector of sun direction (ecef) */
	sunmoonpos(gpst2utc(obs[0].time),erpv,rsun,NULL,NULL);
	normv3(rsun,esun);

	for (i=0;i<n;i++) {
		type=nav->pcvs[obs[i].sat-1].type;

		if ((r=norm(rs+i*6,3))<=0.0) continue;

		/* only block IIA */
		if (*type&&!strstr(type,"BLOCK IIA")) continue;

		/* sun-earth-satellite angle */
		cosa=dot(rs+i*6,esun,3)/r;
		cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);
		ang=acos(cosa);

		/* test eclipse */
		if (ang<PI/2.0||r*sin(ang)>RE_WGS84) continue;

		trace(3,"eclipsing sat excluded %s sat=%2d\n",time_str(obs[0].time,0),
			obs[i].sat);

		for (j=0;j<3;j++) rs[j+i*6]=0.0;
	}
}