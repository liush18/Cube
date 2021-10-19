
#include "base.h"

static char *obscodes[]={       /* observation code strings */

    ""  ,"1C","1P","1W","1Y", "1M","1N","1S","1L","1E", /*  0- 9 */
    "1A","1B","1X","1Z","2C", "2D","2S","2L","2X","2P", /* 10-19 */
    "2W","2Y","2M","2N","5I", "5Q","5X","7I","7Q","7X", /* 20-29 */
    "6A","6B","6C","6X","6Z", "6S","6L","8L","8Q","8X", /* 30-39 */
    "2I","2Q","6I","6Q","3I", "3Q","3X","1I","1Q","5A", /* 40-49 */
    "5B","5C","9A","9B","9C", "9X","1D","5D","5P","5Z", /* 50-59 */
    "6E","7D","7P","7Z","8D", "8P","4A","4B","4X",""    /* 60-69 */
};
static char codepris[7][MAXFREQ][16]={  /* code priority for each freq-index */
    /*    0         1          2          3         4         5     */
    {"CPYWMNSL","PYWCMNDLSX","IQX"     ,""       ,""       ,""      ,""}, /* GPS */
    {"CPABX"   ,"PCABX"     ,"IQX"     ,""       ,""       ,""      ,""}, /* GLO */
    {"CABXZ"   ,"IQX"       ,"IQX"     ,"ABCXZ"  ,"IQX"    ,""      ,""}, /* GAL */
    {"CLSXZ"   ,"LSX"       ,"IQXDPZ"  ,"LSXEZ"  ,""       ,""      ,""}, /* QZS */
    {"C"       ,"IQX"       ,""        ,""       ,""       ,""      ,""}, /* SBS */
    {"IQXDPAN" ,"IQXDPZ"    ,"DPX"     ,"IQXA"   ,"DPX"    ,""      ,""}, /* BDS */
    {"ABCX"    ,"ABCX"      ,""        ,""       ,""       ,""      ,""}  /* IRN */
};


/* test SNR mask ---------------------------------------------------------------
* test SNR mask
* args   : int    base      I   rover or base-station (0:rover,1:base station)
*          int    idx       I   frequency index (0:L1,1:L2,2:L3,...)
*          double el        I   elevation angle (rad)
*          double snr       I   C/N0 (dBHz)
*          snrmask_t *mask  I   SNR mask
* return : status (1:masked,0:unmasked)
*-----------------------------------------------------------------------------*/
extern int testsnr(int base, int idx, double el, double snr, 
    const snrmask_t *mask)
{
    double minsnr,a;
    int i;

    if (!mask->ena[base]||idx<0||idx>=NFREQ) return 0;

    a=(el*R2D+5.0)/10.0;
    i=(int)floor(a); a-=i;
    if      (i<1) minsnr=mask->mask[idx][0];
    else if (i>8) minsnr=mask->mask[idx][8];
    else minsnr=(1.0-a)*mask->mask[idx][i-1]+a*mask->mask[idx][i];

    return snr<minsnr;
}
/* obs type string to obs code -------------------------------------------------
* convert obs code type string to obs code
* args   : char   *str      I   obs code string ("1C","1P","1Y",...)
* return : obs code (CODE_???)
* notes  : obs codes are based on RINEX 3.04
*-----------------------------------------------------------------------------*/
extern uint8_t obs2code(const char *obs)
{
    int i;

    for (i=1;*obscodes[i];i++) {
        if (strcmp(obscodes[i],obs)) continue;
        return (uint8_t)i;
    }
    return CODE_NONE;
}
/* obs code to obs code string -------------------------------------------------
* convert obs code to obs code string
* args   : uint8_t code     I   obs code (CODE_???)
* return : obs code string ("1C","1P","1P",...)
* notes  : obs codes are based on RINEX 3.04
*-----------------------------------------------------------------------------*/
extern char *code2obs(uint8_t code)
{
    if (code<=CODE_NONE||MAXCODE<code) return "";
    return obscodes[code];
}
/* GPS obs code to frequency -------------------------------------------------*/
static int code2freq_GPS(uint8_t code, double *freq)
{
    char *obs=code2obs(code);

    switch (obs[0]) {
        case '1': *freq=FREQ1; return 0; /* L1 */
        case '2': *freq=FREQ2; return 1; /* L2 */
        case '5': *freq=FREQ5; return 2; /* L5 */
    }
    return -1;
}
/* GLONASS obs code to frequency ---------------------------------------------*/
static int code2freq_GLO(uint8_t code, int fcn, double *freq)
{
    char *obs=code2obs(code);

    if (fcn<-7||fcn>6) return -1;

    switch (obs[0]) {
        case '1': *freq=FREQ1_GLO+DFRQ1_GLO*fcn; return 0; /* G1 */
        case '2': *freq=FREQ2_GLO+DFRQ2_GLO*fcn; return 1; /* G2 */
        case '3': *freq=FREQ3_GLO;               return 2; /* G3 */
        case '4': *freq=FREQ1a_GLO;              return 0; /* G1a */
        case '6': *freq=FREQ2a_GLO;              return 1; /* G2a */
    }
    return -1;
}
/* Galileo obs code to frequency ---------------------------------------------*/
static int code2freq_GAL(uint8_t code, double *freq)
{
    char *obs=code2obs(code);

    switch (obs[0]) {
        case '1': *freq=FREQ1; return 0; /* E1 */
        case '7': *freq=FREQ7; return 1; /* E5b */
        case '5': *freq=FREQ5; return 2; /* E5a */
        case '6': *freq=FREQ6; return 3; /* E6 */
        case '8': *freq=FREQ8; return 4; /* E5ab */
    }
    return -1;
}
/* QZSS obs code to frequency ------------------------------------------------*/
static int code2freq_QZS(uint8_t code, double *freq)
{
    char *obs=code2obs(code);

    switch (obs[0]) {
        case '1': *freq=FREQ1; return 0; /* L1 */
        case '2': *freq=FREQ2; return 1; /* L2 */
        case '5': *freq=FREQ5; return 2; /* L5 */
        case '6': *freq=FREQ6; return 3; /* L6 */
    }
    return -1;
}
/* SBAS obs code to frequency ------------------------------------------------*/
static int code2freq_SBS(uint8_t code, double *freq)
{
    char *obs=code2obs(code);

    switch (obs[0]) {
        case '1': *freq=FREQ1; return 0; /* L1 */
        case '5': *freq=FREQ5; return 1; /* L5 */
    }
    return -1;
}
/* BDS obs code to frequency -------------------------------------------------*/
static int code2freq_BDS(uint8_t code, double *freq)
{
    char *obs=code2obs(code);

    switch (obs[0]) {
        case '1': *freq=FREQ1;     return 0; /* B1C */
        case '2': *freq=FREQ1_CMP; return 0; /* B1I */
        case '7': *freq=FREQ2_CMP; return 1; /* B2I/B2b */
        case '5': *freq=FREQ5;     return 2; /* B2a */
        case '6': *freq=FREQ3_CMP; return 3; /* B3 */
        case '8': *freq=FREQ8;     return 4; /* B2ab */
    }
    return -1;
}
/* NavIC obs code to frequency -----------------------------------------------*/
static int code2freq_IRN(uint8_t code, double *freq)
{
    char *obs=code2obs(code);

    switch (obs[0]) {
        case '5': *freq=FREQ5; return 0; /* L5 */
        case '9': *freq=FREQ9; return 1; /* S */
    }
    return -1;
}
/* system and obs code to frequency index --------------------------------------
* convert system and obs code to frequency index
* args   : int    sys       I   satellite system (SYS_???)
*          uint8_t code     I   obs code (CODE_???)
* return : frequency index (-1: error)
*                       0     1     2     3     4 
*           --------------------------------------
*            GPS       L1    L2    L5     -     - 
*            GLONASS   G1    G2    G3     -     -  (G1=G1,G1a,G2=G2,G2a)
*            Galileo   E1    E5b   E5a   E6   E5ab
*            QZSS      L1    L2    L5    L6     - 
*            SBAS      L1     -    L5     -     -
*            BDS       B1    B2    B2a   B3   B2ab (B1=B1I,B1C,B2=B2I,B2b)
*            NavIC     L5     S     -     -     - 
*-----------------------------------------------------------------------------*/
extern int code2idx(int sys, uint8_t code)
{
    double freq;

    switch (sys) {
        case SYS_GPS: return code2freq_GPS(code,&freq);
        case SYS_GLO: return code2freq_GLO(code,0,&freq);
        case SYS_GAL: return code2freq_GAL(code,&freq);
        case SYS_QZS: return code2freq_QZS(code,&freq);
        case SYS_SBS: return code2freq_SBS(code,&freq);
        case SYS_CMP: return code2freq_BDS(code,&freq);
        case SYS_IRN: return code2freq_IRN(code,&freq);
    }
    return -1;
}
/* system and obs code to frequency --------------------------------------------
* convert system and obs code to carrier frequency
* args   : int    sys       I   satellite system (SYS_???)
*          uint8_t code     I   obs code (CODE_???)
*          int    fcn       I   frequency channel number for GLONASS
* return : carrier frequency (Hz) (0.0: error)
*-----------------------------------------------------------------------------*/
extern double code2freq(int sys, uint8_t code, int fcn)
{
    double freq=0.0;

    switch (sys) {
        case SYS_GPS: (void)code2freq_GPS(code,&freq); break;
        case SYS_GLO: (void)code2freq_GLO(code,fcn,&freq); break;
        case SYS_GAL: (void)code2freq_GAL(code,&freq); break;
        case SYS_QZS: (void)code2freq_QZS(code,&freq); break;
        case SYS_SBS: (void)code2freq_SBS(code,&freq); break;
        case SYS_CMP: (void)code2freq_BDS(code,&freq); break;
        case SYS_IRN: (void)code2freq_IRN(code,&freq); break;
    }
    return freq;
}
/* satellite and obs code to frequency -----------------------------------------
* convert satellite and obs code to carrier frequency
* args   : int    sat       I   satellite number
*          uint8_t code     I   obs code (CODE_???)
*          nav_t  *nav_t    I   navigation data for GLONASS (NULL: not used)
* return : carrier frequency (Hz) (0.0: error)
*-----------------------------------------------------------------------------*/
extern double sat2freq(int sat, uint8_t code, const nav_t *nav)
{
    int i,fcn=0,sys,prn;

    sys=satsys(sat,&prn);

    if (sys==SYS_GLO) {
        if (!nav) return 0.0;
        for (i=0;i<nav->ng;i++) {
            if (nav->geph[i].sat==sat) break;
        }
        if (i<nav->ng) {
            fcn=nav->geph[i].frq;
        }
        else if (nav->glo_fcn[prn-1]>0) {
            fcn=nav->glo_fcn[prn-1]-8;
        }
        else return 0.0;
    }
    return code2freq(sys,code,fcn);
}
/* set code priority -----------------------------------------------------------
* set code priority for multiple codes in a frequency
* args   : int    sys       I   system (or of SYS_???)
*          int    idx       I   frequency index (0- )
*          char   *pri      I   priority of codes (series of code characters)
*                               (higher priority precedes lower)
* return : none
*-----------------------------------------------------------------------------*/
extern void setcodepri(int sys, int idx, const char *pri)
{
    trace(3,"setcodepri:sys=%d idx=%d pri=%s\n",sys,idx,pri);

    if (idx<0||idx>=MAXFREQ) return;
    if (sys&SYS_GPS) strcpy(codepris[0][idx],pri);
    if (sys&SYS_GLO) strcpy(codepris[1][idx],pri);
    if (sys&SYS_GAL) strcpy(codepris[2][idx],pri);
    if (sys&SYS_QZS) strcpy(codepris[3][idx],pri);
    if (sys&SYS_SBS) strcpy(codepris[4][idx],pri);
    if (sys&SYS_CMP) strcpy(codepris[5][idx],pri);
    if (sys&SYS_IRN) strcpy(codepris[6][idx],pri);
}
/* get code priority -----------------------------------------------------------
* get code priority for multiple codes in a frequency
* args   : int    sys       I   system (SYS_???)
*          uint8_t code     I   obs code (CODE_???)
*          char   *opt      I   code options (NULL:no option)
* return : priority (15:highest-1:lowest,0:error)
*-----------------------------------------------------------------------------*/
extern int getcodepri(int sys, uint8_t code, const char *opt)
{
    const char *p,*optstr;
    char *obs,str[8]="";
    int i,j;

    switch (sys) {
        case SYS_GPS: i=0; optstr="-GL%2s"; break;
        case SYS_GLO: i=1; optstr="-RL%2s"; break;
        case SYS_GAL: i=2; optstr="-EL%2s"; break;
        case SYS_QZS: i=3; optstr="-JL%2s"; break;
        case SYS_SBS: i=4; optstr="-SL%2s"; break;
        case SYS_CMP: i=5; optstr="-CL%2s"; break;
        case SYS_IRN: i=6; optstr="-IL%2s"; break;
        default: return 0;
    }
    if ((j=code2idx(sys,code))<0) return 0;
    obs=code2obs(code);

    /* parse code options */
    for (p=opt;p&&(p=strchr(p,'-'));p++) {
        if (sscanf(p,optstr,str)<1||str[0]!=obs[0]) continue;
        return str[1]==obs[1]?15:0;
    }
    /* search code priority */
    return (p=strchr(codepris[i][j],obs[1]))?14-(int)(p-codepris[i][j]):0;
}
/* carrier-phase LC (m) ------------------------------------------------------*/
static double L_LC(int i,int j,int k,const double *L)
{
	const double f1=FREQ1,f2=FREQ2,f5=FREQ5;
	double L1,L2,L5;
	if ((i&&!L[0])||(j&&!L[1])||(k&&!L[2])) return 0.0;
	L1=L[0];
	L2=L[1];
	L5=L[2];
	return (i*f1*L1+j*f2*L2+k*f5*L5)/(i*f1+j*f2+k*f5);
}
/* pseudorange LC (m) --------------------------------------------------------*/
static double P_LC(int i,int j,int k,const double *P)
{
	const double f1=FREQ1,f2=FREQ2,f5=FREQ5;
	double P1,P2,P5;

	if ((i&&!P[0])||(j&&!P[1])||(k&&!P[2])) return 0.0;
	P1=P[0];
	P2=P[1];
	P5=P[2];
	return (i*f1*P1+j*f2*P2+k*f5*P5)/(i*f1+j*f2+k*f5);
}
/* Melbourne-W¨¹bbena widelane combination ------------------------------------*/
extern void combWideLane(const double *L,const double *P,double *A4)
{
	*A4=L_LC(1,-1,0,L)-P_LC(1,1,0,P);
}
#if 0
/* Melbourne-W¨¹bbena widelane combination ------------------------------------*/
extern void calWideLane(const obsd_t *obs, const nav_t *nav, const double *azel,
	const prcopt_t *opt, double *A4)
{
	const double *lam=nav->lam[obs->sat-1];
	double L[NFREQ],P[NFREQ];
	int i;

	for (i=0;i<NFREQ;i++) {
		L[i]=P[i]=0.0;
		if (lam[i]==0.0||obs->L[i]==0.0||obs->P[i]==0.0) continue;
		if (testsnr(0,0,azel[1],obs->SNR[i]*0.25,&opt->snrmask)) continue;

		/* antenna phase center and phase windup correction */
		L[i]=obs->L[i]*lam[i];
		P[i]=obs->P[i];

		/* P1-C1,P2-C2 dcb correction (C1->P1,C2->P2) */
		if (obs->code[i]==CODE_L1C) {
			P[i]+=nav->cbias[obs->sat-1][1];
		}
		else if (obs->code[i]==CODE_L2C||obs->code[i]==CODE_L2X||
			obs->code[i]==CODE_L2L||obs->code[i]==CODE_L2S) {
			P[i]+=nav->cbias[obs->sat-1][2];
		}
	}
	*A4=0.0;
	if (L[0]!=0.0&&L[1]!=0.0&&P[0]!=0.0&&P[1]!=0.0) 
		*A4=L_LC(1,-1,0,L)-P_LC(1,1,0,P);
}
#endif