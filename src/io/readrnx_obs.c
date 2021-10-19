
#include "io.h"


static const char syscodes[]="GREJSCI"; /* satellite system codes */
static const char obscodes[]="CLDS";    /* observation type codes */

/* type definition -----------------------------------------------------------*/
typedef struct {                        /* signal index type */
    int n;                              /* number of index */
    int idx[MAXOBSTYPE];                /* signal freq-index */
    int pos[MAXOBSTYPE];                /* signal index in obs data (-1:no) */
    uint8_t pri [MAXOBSTYPE];           /* signal priority (15-0) */
    uint8_t type[MAXOBSTYPE];           /* type (0:C,1:L,2:D,3:S) */
    uint8_t code[MAXOBSTYPE];           /* obs-code (CODE_L??) */
    double shift[MAXOBSTYPE];           /* phase shift (cycle) */
} sigind_t;

/* convert RINEX obs-type ver.2 -> ver.3 -------------------------------------*/
static void convcode(double ver, int sys, const char *str, char *type)
{
    strcpy(type,"   ");

    if      (!strcmp(str,"P1")) { /* ver.2.11 GPS L1PY,GLO L2P */
        if      (sys==SYS_GPS) sprintf(type,"%c1W",'C');
        else if (sys==SYS_GLO) sprintf(type,"%c1P",'C');
    }
    else if (!strcmp(str,"P2")) { /* ver.2.11 GPS L2PY,GLO L2P */
        if      (sys==SYS_GPS) sprintf(type,"%c2W",'C');
        else if (sys==SYS_GLO) sprintf(type,"%c2P",'C');
    }
    else if (!strcmp(str,"C1")) { /* ver.2.11 GPS L1C,GLO L1C/A */
        if      (ver>=2.12) ; /* reject C1 for 2.12 */
        else if (sys==SYS_GPS) sprintf(type,"%c1C",'C');
        else if (sys==SYS_GLO) sprintf(type,"%c1C",'C');
        else if (sys==SYS_GAL) sprintf(type,"%c1X",'C'); /* ver.2.12 */
        else if (sys==SYS_QZS) sprintf(type,"%c1C",'C');
        else if (sys==SYS_SBS) sprintf(type,"%c1C",'C');
    }
    else if (!strcmp(str,"C2")) {
        if (sys==SYS_GPS) {
            if (ver>=2.12) sprintf(type,"%c2W",'C'); /* L2P(Y) */
            else           sprintf(type,"%c2X",'C'); /* L2C */
        }
        else if (sys==SYS_GLO) sprintf(type,"%c2C",'C');
        else if (sys==SYS_QZS) sprintf(type,"%c2X",'C');
        else if (sys==SYS_CMP) sprintf(type,"%c2X",'C'); /* ver.2.12 B1_2 */
    }
    else if (ver>=2.12&&str[1]=='A') { /* ver.2.12 L1C/A */
        if      (sys==SYS_GPS) sprintf(type,"%c1C",str[0]);
        else if (sys==SYS_GLO) sprintf(type,"%c1C",str[0]);
        else if (sys==SYS_QZS) sprintf(type,"%c1C",str[0]);
        else if (sys==SYS_SBS) sprintf(type,"%c1C",str[0]);
    }
    else if (ver>=2.12&&str[1]=='B') { /* ver.2.12 GPS L1C */
        if      (sys==SYS_GPS) sprintf(type,"%c1X",str[0]);
        else if (sys==SYS_QZS) sprintf(type,"%c1X",str[0]);
    }
    else if (ver>=2.12&&str[1]=='C') { /* ver.2.12 GPS L2C */
        if      (sys==SYS_GPS) sprintf(type,"%c2X",str[0]);
        else if (sys==SYS_QZS) sprintf(type,"%c2X",str[0]);
    }
    else if (ver>=2.12&&str[1]=='D') { /* ver.2.12 GLO L2C/A */
        if      (sys==SYS_GLO) sprintf(type,"%c2C",str[0]);
    }
    else if (ver>=2.12&&str[1]=='1') { /* ver.2.12 GPS L1PY,GLO L1P */
        if      (sys==SYS_GPS) sprintf(type,"%c1W",str[0]);
        else if (sys==SYS_GLO) sprintf(type,"%c1P",str[0]);
        else if (sys==SYS_GAL) sprintf(type,"%c1X",str[0]); /* tentative */
        else if (sys==SYS_CMP) sprintf(type,"%c2X",str[0]); /* extension */
    }
    else if (ver<2.12&&str[1]=='1') {
        if      (sys==SYS_GPS) sprintf(type,"%c1C",str[0]);
        else if (sys==SYS_GLO) sprintf(type,"%c1C",str[0]);
        else if (sys==SYS_GAL) sprintf(type,"%c1X",str[0]); /* tentative */
        else if (sys==SYS_QZS) sprintf(type,"%c1C",str[0]);
        else if (sys==SYS_SBS) sprintf(type,"%c1C",str[0]);
    }
    else if (str[1]=='2') {
        if      (sys==SYS_GPS) sprintf(type,"%c2W",str[0]);
        else if (sys==SYS_GLO) sprintf(type,"%c2P",str[0]);
        else if (sys==SYS_QZS) sprintf(type,"%c2X",str[0]);
        else if (sys==SYS_CMP) sprintf(type,"%c2X",str[0]); /* ver.2.12 B1_2 */
    }
    else if (str[1]=='5') {
        if      (sys==SYS_GPS) sprintf(type,"%c5X",str[0]);
        else if (sys==SYS_GAL) sprintf(type,"%c5X",str[0]);
        else if (sys==SYS_QZS) sprintf(type,"%c5X",str[0]);
        else if (sys==SYS_SBS) sprintf(type,"%c5X",str[0]);
    }
    else if (str[1]=='6') {
        if      (sys==SYS_GAL) sprintf(type,"%c6X",str[0]);
        else if (sys==SYS_QZS) sprintf(type,"%c6X",str[0]);
        else if (sys==SYS_CMP) sprintf(type,"%c6X",str[0]); /* ver.2.12 B3 */
    }
    else if (str[1]=='7') {
        if      (sys==SYS_GAL) sprintf(type,"%c7X",str[0]);
        else if (sys==SYS_CMP) sprintf(type,"%c7X",str[0]); /* ver.2.12 B2b */
    }
    else if (str[1]=='8') {
        if      (sys==SYS_GAL) sprintf(type,"%c8X",str[0]);
    }
    trace(3,"convcode: ver=%.2f sys=%2d type= %s -> %s\n",ver,sys,str,type);
}
/* decode RINEX observation data file header ---------------------------------*/
extern void decode_obsh(FILE *fp, char *buff, double ver, int *tsys,
    char tobs[][MAXOBSTYPE][4], nav_t *nav, sta_t *sta)
{
    /* default codes for unknown code */
    const char frqcodes[]="1256789";
    const char *defcodes[]={
        "CWX    ",  /* GPS: L125____ */
        "CCXX X ",  /* GLO: L1234_6_ */
        "C XXXX ",  /* GAL: L1_5678_ */
        "CXXX   ",  /* QZS: L1256___ */
        "C X    ",  /* SBS: L1_5____ */
        "XIXIIX ",  /* BDS: L125678_ */
        "  A   A"   /* IRN: L__5___9 */
    };
    double del[3];
    int i,j,k,n,nt,prn,fcn;
    const char *p;
    char *label=buff+60,str[4];

    trace(4,"decode_obsh: ver=%.2f\n",ver);

    if      (strstr(label,"MARKER NAME"         )) {
        if (sta) setstr(sta->name,buff,60);
    }
    else if (strstr(label,"MARKER NUMBER"       )) { /* opt */
        if (sta) setstr(sta->marker,buff,20);
    }
    else if (strstr(label,"MARKER TYPE"         )) ; /* ver.3 */
    else if (strstr(label,"OBSERVER / AGENCY"   )) ;
    else if (strstr(label,"REC # / TYPE / VERS" )) {
        if (sta) {
            setstr(sta->recsno, buff,   20);
            setstr(sta->rectype,buff+20,20);
            setstr(sta->recver, buff+40,20);
        }
    }
    else if (strstr(label,"ANT # / TYPE"        )) {
        if (sta) {
            setstr(sta->antsno,buff   ,20);
            setstr(sta->antdes,buff+20,20);
        }
    }
    else if (strstr(label,"APPROX POSITION XYZ" )) {
        if (sta) {
            for (i=0,j=0;i<3;i++,j+=14) sta->pos[i]=str2num(buff,j,14);
        }
    }
    else if (strstr(label,"ANTENNA: DELTA H/E/N")) {
        if (sta) {
            for (i=0,j=0;i<3;i++,j+=14) del[i]=str2num(buff,j,14);
            sta->del[2]=del[0]; /* h */
            sta->del[0]=del[1]; /* e */
            sta->del[1]=del[2]; /* n */
        }
    }
    else if (strstr(label,"ANTENNA: DELTA X/Y/Z")) ; /* opt ver.3 */
    else if (strstr(label,"ANTENNA: PHASECENTER")) ; /* opt ver.3 */
    else if (strstr(label,"ANTENNA: B.SIGHT XYZ")) ; /* opt ver.3 */
    else if (strstr(label,"ANTENNA: ZERODIR AZI")) ; /* opt ver.3 */
    else if (strstr(label,"ANTENNA: ZERODIR XYZ")) ; /* opt ver.3 */
    else if (strstr(label,"CENTER OF MASS: XYZ" )) ; /* opt ver.3 */
    else if (strstr(label,"SYS / # / OBS TYPES" )) { /* ver.3 */
        if (!(p=strchr(syscodes,buff[0]))) {
            trace(2,"invalid system code: sys=%c\n",buff[0]);
            return;
        }
        i=(int)(p-syscodes);
        n=(int)str2num(buff,3,3);
        for (j=nt=0,k=7;j<n;j++,k+=4) {
            if (k>58) {
                if (!fgets(buff,MAXRNXLEN,fp)) break;
                k=7;
            }
            if (nt<MAXOBSTYPE-1) setstr(tobs[i][nt++],buff+k,3);
        }
        *tobs[i][nt]='\0';

        /* change BDS B1 code: 3.02 */
        if (i==5&&fabs(ver-3.02)<1e-3) {
            for (j=0;j<nt;j++) if (tobs[i][j][1]=='1') tobs[i][j][1]='2';
        }
        /* if unknown code in ver.3, set default code */
        for (j=0;j<nt;j++) {
            if (tobs[i][j][2]) continue;
            if (!(p=strchr(frqcodes,tobs[i][j][1]))) continue;
            tobs[i][j][2]=defcodes[i][(int)(p-frqcodes)];
            trace(2,"set default for unknown code: sys=%c code=%s\n",buff[0],
                tobs[i][j]);
        }
    }
    else if (strstr(label,"WAVELENGTH FACT L1/2")) ; /* opt ver.2 */
    else if (strstr(label,"# / TYPES OF OBSERV" )) { /* ver.2 */
        n=(int)str2num(buff,0,6);
        for (i=nt=0,j=10;i<n;i++,j+=6) {
            if (j>58) {
                if (!fgets(buff,MAXRNXLEN,fp)) break;
                j=10;
            }
            if (nt>=MAXOBSTYPE-1) continue;
            if (ver<=2.99) {
                setstr(str,buff+j,2);
                convcode(ver,SYS_GPS,str,tobs[0][nt]);
                convcode(ver,SYS_GLO,str,tobs[1][nt]);
                convcode(ver,SYS_GAL,str,tobs[2][nt]);
                convcode(ver,SYS_QZS,str,tobs[3][nt]);
                convcode(ver,SYS_SBS,str,tobs[4][nt]);
                convcode(ver,SYS_CMP,str,tobs[5][nt]);
            }
            nt++;
        }
        *tobs[0][nt]='\0';
    }
    else if (strstr(label,"SIGNAL STRENGTH UNIT")) ; /* opt ver.3 */
    else if (strstr(label,"INTERVAL"            )) ; /* opt */
    else if (strstr(label,"TIME OF FIRST OBS"   )) {
        if      (!strncmp(buff+48,"GPS",3)) *tsys=TSYS_GPS;
        else if (!strncmp(buff+48,"GLO",3)) *tsys=TSYS_UTC;
        else if (!strncmp(buff+48,"GAL",3)) *tsys=TSYS_GAL;
        else if (!strncmp(buff+48,"QZS",3)) *tsys=TSYS_QZS; /* ver.3.02 */
        else if (!strncmp(buff+48,"BDT",3)) *tsys=TSYS_CMP; /* ver.3.02 */
        else if (!strncmp(buff+48,"IRN",3)) *tsys=TSYS_IRN; /* ver.3.03 */
    }
    else if (strstr(label,"TIME OF LAST OBS"    )) ; /* opt */
    else if (strstr(label,"RCV CLOCK OFFS APPL" )) ; /* opt */
    else if (strstr(label,"SYS / DCBS APPLIED"  )) ; /* opt ver.3 */
    else if (strstr(label,"SYS / PCVS APPLIED"  )) ; /* opt ver.3 */
    else if (strstr(label,"SYS / SCALE FACTOR"  )) ; /* opt ver.3 */
    else if (strstr(label,"SYS / PHASE SHIFTS"  )) ; /* ver.3.01 */
    else if (strstr(label,"GLONASS SLOT / FRQ #")) { /* ver.3.02 */
        for (i=0;i<8;i++) {
            if (buff[4+i*7]!='R') continue;
            prn=(int)str2num(buff,5+i*7,2);
            fcn=(int)str2num(buff,8+i*7,2);
            if (prn<1||prn>MAXPRNGLO||fcn<-7||fcn>6) continue;
            if (nav) nav->glo_fcn[prn-1]=fcn+8;
        }
    }
    else if (strstr(label,"GLONASS COD/PHS/BIS" )) { /* ver.3.02 */
        if (sta) {
            sta->glo_cp_bias[0]=str2num(buff, 5,8);
            sta->glo_cp_bias[1]=str2num(buff,18,8);
            sta->glo_cp_bias[2]=str2num(buff,31,8);
            sta->glo_cp_bias[3]=str2num(buff,44,8);
        }
    }
    else if (strstr(label,"LEAP SECONDS"        )) { /* opt */
        if (nav) {
            nav->utc_gps[4]=str2num(buff, 0,6);
            nav->utc_gps[7]=str2num(buff, 6,6);
            nav->utc_gps[5]=str2num(buff,12,6);
            nav->utc_gps[6]=str2num(buff,18,6);
        }
    }
    else if (strstr(label,"# OF SALTELLITES"    )) { /* opt */
        /* skip */ ;
    }
    else if (strstr(label,"PRN / # OF OBS"      )) { /* opt */
        /* skip */ ;
    }
}
/* decode observation epoch --------------------------------------------------*/
static int decode_obsepoch(FILE *fp, char *buff, double ver, gtime_t *time,
    int *flag, int *sats)
{
    int i,j,n;
    char satid[8]="";

    trace(4,"decode_obsepoch: ver=%.2f\n",ver);

    if (ver<=2.99) { /* ver.2 */
        if ((n=(int)str2num(buff,29,3))<=0) return 0;

        /* epoch flag: 3:new site,4:header info,5:external event */
        *flag=(int)str2num(buff,28,1);

        if (3<=*flag&&*flag<=5) return n;

        if (str2time(buff,0,26,time)) {
            trace(2,"rinex obs invalid epoch: epoch=%26.26s\n",buff);
            return 0;
        }
        for (i=0,j=32;i<n;i++,j+=3) {
            if (j>=68) {
                if (!fgets(buff,MAXRNXLEN,fp)) break;
                j=32;
            }
            if (i<MAXOBS) {
                strncpy(satid,buff+j,3);
                sats[i]=satid2no(satid);
            }
        }
    }
    else { /* ver.3 */
        if ((n=(int)str2num(buff,32,3))<=0) return 0;

        *flag=(int)str2num(buff,31,1);

        if (3<=*flag&&*flag<=5) return n;

        if (buff[0]!='>'||str2time(buff,1,28,time)) {
            trace(2,"rinex obs invalid epoch: epoch=%29.29s\n",buff);
            return 0;
        }
    }
    trace(4,"decode_obsepoch: time=%s flag=%d\n",time_str(*time,3),*flag);
    return n;
}
/* decode observation data ---------------------------------------------------*/
static int decode_obsdata(FILE *fp, char *buff, double ver, int mask,
    sigind_t *index, obsd_t *obs)
{
    sigind_t *ind;
    double val[MAXOBSTYPE]={0};
    uint8_t lli[MAXOBSTYPE]={0};
    char satid[8]="";
    int i,j,n,m,stat=1,p[MAXOBSTYPE],k[16],l[16];

    trace(4,"decode_obsdata: ver=%.2f\n",ver);

    if (ver>2.99) { /* ver.3 */
        sprintf(satid,"%.3s",buff);
        obs->sat=(uint8_t)satid2no(satid);
    }
    if (!obs->sat) {
        trace(4,"decode_obsdata: unsupported sat sat=%s\n",satid);
        stat=0;
    }
    else if (!(satsys(obs->sat,NULL)&mask)) {
        stat=0;
    }
    /* read observation data fields */
    switch (satsys(obs->sat,NULL)) {
        case SYS_GLO: ind=index+1; break;
        case SYS_GAL: ind=index+2; break;
        case SYS_QZS: ind=index+3; break;
        case SYS_SBS: ind=index+4; break;
        case SYS_CMP: ind=index+5; break;
        case SYS_IRN: ind=index+6; break;
        default:      ind=index  ; break;
    }
    for (i=0,j=ver<=2.99?0:3;i<ind->n;i++,j+=16) {

        if (ver<=2.99&&j>=80) { /* ver.2 */
            if (!fgets(buff,MAXRNXLEN,fp)) break;
            j=0;
        }
        if (stat) {
            val[i]=str2num(buff,j,14)+ind->shift[i];
            lli[i]=(uint8_t)str2num(buff,j+14,1)&3;
        }
    }
    if (!stat) return 0;

    for (i=0;i<NFREQ+NEXOBS;i++) {
        obs->P[i]=obs->L[i]=0.0; obs->D[i]=0.0f;
        obs->SNR[i]=obs->LLI[i]=obs->code[i]=0;
    }
    /* assign position in observation data */
    for (i=n=m=0;i<ind->n;i++) {

        p[i]=(ver<=2.11)?ind->idx[i]:ind->pos[i];

        if (ind->type[i]==0&&p[i]==0) k[n++]=i; /* C1? index */
        if (ind->type[i]==0&&p[i]==1) l[m++]=i; /* C2? index */
    }
    if (ver<=2.11) {

        /* if multiple codes (C1/P1,C2/P2), select higher priority */
        if (n>=2) {
            if (val[k[0]]==0.0&&val[k[1]]==0.0) {
                p[k[0]]=-1; p[k[1]]=-1;
            }
            else if (val[k[0]]!=0.0&&val[k[1]]==0.0) {
                p[k[0]]=0; p[k[1]]=-1;
            }
            else if (val[k[0]]==0.0&&val[k[1]]!=0.0) {
                p[k[0]]=-1; p[k[1]]=0;
            }
            else if (ind->pri[k[1]]>ind->pri[k[0]]) {
                p[k[1]]=0; p[k[0]]=NEXOBS<1?-1:NFREQ;
            }
            else {
                p[k[0]]=0; p[k[1]]=NEXOBS<1?-1:NFREQ;
            }
        }
        if (m>=2) {
            if (val[l[0]]==0.0&&val[l[1]]==0.0) {
                p[l[0]]=-1; p[l[1]]=-1;
            }
            else if (val[l[0]]!=0.0&&val[l[1]]==0.0) {
                p[l[0]]=1; p[l[1]]=-1;
            }
            else if (val[l[0]]==0.0&&val[l[1]]!=0.0) {
                p[l[0]]=-1; p[l[1]]=1; 
            }
            else if (ind->pri[l[1]]>ind->pri[l[0]]) {
                p[l[1]]=1; p[l[0]]=NEXOBS<2?-1:NFREQ+1;
            }
            else {
                p[l[0]]=1; p[l[1]]=NEXOBS<2?-1:NFREQ+1;
            }
        }
    }
    /* save observation data */
    for (i=0;i<ind->n;i++) {
        if (p[i]<0||val[i]==0.0) continue;
        switch (ind->type[i]) {
        case 0: obs->P[p[i]]=val[i]; obs->code[p[i]]=ind->code[i]; break;
        case 1: obs->L[p[i]]=val[i]; obs->LLI [p[i]]=lli[i];    break;
        case 2: obs->D[p[i]]=(float)val[i];                     break;
        case 3: obs->SNR[p[i]]=(uint16_t)(val[i]/SNR_UNIT+0.5); break;
        }
    }
    trace(4,"decode_obsdata: time=%s sat=%2d\n",time_str(obs->time,0),obs->sat);
    return 1;
}
/* save cycle slips ----------------------------------------------------------*/
extern void saveslips(uint8_t slips[][NFREQ+NEXOBS], obsd_t *data)
{
    int i;
    for (i=0;i<NFREQ+NEXOBS;i++) {
        if (data->LLI[i]&1) slips[data->sat-1][i]|=LLI_SLIP;
    }
}
/* restore cycle slips -------------------------------------------------------*/
extern void restslips(uint8_t slips[][NFREQ+NEXOBS], obsd_t *data)
{
    int i;
    for (i=0;i<NFREQ+NEXOBS;i++) {
        if (slips[data->sat-1][i]&1) data->LLI[i]|=LLI_SLIP;
        slips[data->sat-1][i]=0;
    }
}
/* add observation data ------------------------------------------------------*/
static int addobsdata(obs_t *obs, const obsd_t *data)
{
    obsd_t *obs_data;

    if (obs->nmax<=obs->n) {
        if (obs->nmax<=0) obs->nmax=NINCOBS; else obs->nmax*=2;
        if (!(obs_data=(obsd_t *)realloc(obs->data,sizeof(obsd_t)*obs->nmax))) {
            trace(1,"addobsdata: malloc error n=%dx%d\n",sizeof(obsd_t),obs->nmax);
            free(obs->data); obs->data=NULL; obs->n=obs->nmax=0;
            return -1;
        }
        obs->data=obs_data;
    }
    obs->data[obs->n++]=*data;
    return 1;
}
/* set signal index ----------------------------------------------------------*/
static void set_index(double ver, int sys, const char *opt,
    char tobs[MAXOBSTYPE][4], sigind_t *ind)
{
    const char *p;
    char str[8],*optstr="";
    double shift;
    int i,j,k,n;

    for (i=n=0;*tobs[i];i++,n++) {
        ind->code[i]=obs2code(tobs[i]+1);
        ind->type[i]=(p=strchr(obscodes,tobs[i][0]))?(int)(p-obscodes):0;
        ind->idx[i]=code2idx(sys,ind->code[i]);
        ind->pri[i]=getcodepri(sys,ind->code[i],opt);
        ind->pos[i]=-1;
    }
    /* parse phase shift options */
    switch (sys) {
        case SYS_GPS: optstr="-GL%2s=%lf"; break;
        case SYS_GLO: optstr="-RL%2s=%lf"; break;
        case SYS_GAL: optstr="-EL%2s=%lf"; break;
        case SYS_QZS: optstr="-JL%2s=%lf"; break;
        case SYS_SBS: optstr="-SL%2s=%lf"; break;
        case SYS_CMP: optstr="-CL%2s=%lf"; break;
        case SYS_IRN: optstr="-IL%2s=%lf"; break;
    }
    for (p=opt;p&&(p=strchr(p,'-'));p++) {
        if (sscanf(p,optstr,str,&shift)<2) continue;
        for (i=0;i<n;i++) {
            if (strcmp(code2obs(ind->code[i]),str)) continue;
            ind->shift[i]=shift;
            trace(2,"phase shift: sys=%2d tobs=%s shift=%.3f\n",sys,
                tobs[i],shift);
        }
    }
    /* assign index for highest priority code */
    for (i=0;i<NFREQ;i++) {
        for (j=0,k=-1;j<n;j++) {
            if (ind->idx[j]==i&&ind->pri[j]&&(k<0||ind->pri[j]>ind->pri[k])) {
                k=j;
            }
        }
        if (k<0) continue;

        for (j=0;j<n;j++) {
            if (ind->code[j]==ind->code[k]) ind->pos[j]=i;
        }
    }
    /* assign index of extended observation data */
    for (i=0;i<NEXOBS;i++) {
        for (j=0;j<n;j++) {
            if (ind->code[j]&&ind->pri[j]&&ind->pos[j]<0) break;
        }
        if (j>=n) break;

        for (k=0;k<n;k++) {
            if (ind->code[k]==ind->code[j]) ind->pos[k]=NFREQ+i;
        }
    }
    for (i=0;i<n;i++) {
        if (!ind->code[i]||!ind->pri[i]||ind->pos[i]>=0) continue;
        trace(4,"reject obs type: sys=%2d, obs=%s\n",sys,tobs[i]);
    }
    ind->n=n;

#if 0 /* for debug */
    for (i=0;i<n;i++) {
        trace(2,"set_index: sys=%2d,tobs=%s code=%2d pri=%2d idx=%d pos=%d shift=%5.2f\n",
            sys,tobs[i],ind->code[i],ind->pri[i],ind->idx[i],ind->pos[i],
            ind->shift[i]);
    }
#endif
}
/* read RINEX observation data body ------------------------------------------*/
extern int readrnxobsb(FILE *fp, const char *opt, double ver, int *tsys,
    char tobs[][MAXOBSTYPE][4], int *flag, obsd_t *data, sta_t *sta)
{
    gtime_t time={0};
    sigind_t index[NUMSYS]={{0}};
    char buff[MAXRNXLEN];
    int i=0,n=0,nsat=0,nsys=NUMSYS,sats[MAXOBS]={0},mask;

    /* set system mask */
    mask=set_sysmask(opt);

    /* set signal index */
    if (nsys>=1) set_index(ver,SYS_GPS,opt,tobs[0],index  );
    if (nsys>=2) set_index(ver,SYS_GLO,opt,tobs[1],index+1);
    if (nsys>=3) set_index(ver,SYS_GAL,opt,tobs[2],index+2);
    if (nsys>=4) set_index(ver,SYS_QZS,opt,tobs[3],index+3);
    if (nsys>=5) set_index(ver,SYS_SBS,opt,tobs[4],index+4);
    if (nsys>=6) set_index(ver,SYS_CMP,opt,tobs[5],index+5);
    if (nsys>=7) set_index(ver,SYS_IRN,opt,tobs[6],index+6);

    /* read record */
    while (fgets(buff,MAXRNXLEN,fp)) {

        /* decode observation epoch */
        if (i==0) {
            if ((nsat=decode_obsepoch(fp,buff,ver,&time,flag,sats))<=0) {
                continue;
            }
        }
        else if ((*flag<=2||*flag==6)&&n<MAXOBS) {
            data[n].time=time;
            data[n].sat=(uint8_t)sats[i-1];

            /* decode RINEX observation data */
            if (decode_obsdata(fp,buff,ver,mask,index,data+n)) n++;
        }
        else if (*flag==3||*flag==4) { /* new site or header info follows */

                                       /* decode RINEX observation data file header */
            decode_obsh(fp,buff,ver,tsys,tobs,NULL,sta);
        }
        if (++i>nsat) return n;
    }
    return -1;
}
/* read RINEX observation data -----------------------------------------------*/
extern int readrnxobs(FILE *fp, gtime_t ts, gtime_t te, double tint,
    const char *opt, int rcv, double ver, int *tsys, char tobs[][MAXOBSTYPE][4],
    obs_t *obs, sta_t *sta)
{
    obsd_t *data;
    uint8_t slips[MAXSAT][NFREQ+NEXOBS]={{0}};
    int i,n,flag=0,stat=0;

    trace(4,"readrnxobs: rcv=%d ver=%.2f tsys=%d\n",rcv,ver,tsys);

    if (!obs||rcv>MAXRCV) return 0;

    if (!(data=(obsd_t *)malloc(sizeof(obsd_t)*MAXOBS))) return 0;

    /* read RINEX observation data body */
    while ((n=readrnxobsb(fp,opt,ver,tsys,tobs,&flag,data,sta))>=0&&stat>=0) {

        for (i=0;i<n;i++) {

            /* UTC -> GPST */
            if (*tsys==TSYS_UTC) data[i].time=utc2gpst(data[i].time);

            /* save cycle slip */
            saveslips(slips,data+i);
        }
        /* screen data by time */
        if (n>0&&!screent(data[0].time,ts,te,tint)) continue;

        for (i=0;i<n;i++) {

            /* restore cycle slip */
            restslips(slips,data+i);

            data[i].rcv=(uint8_t)rcv;

            /* save obs data */
            if ((stat=addobsdata(obs,data+i))<0) break;
        }
    }
    trace(4,"readrnxobs: nobs=%d stat=%d\n",obs->n,stat);

    free(data);

    return stat;
}