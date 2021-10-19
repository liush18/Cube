
#include "io.h"

/* decode antenna parameter field --------------------------------------------*/
static int decodef(char* p, int n, double* v)
{
    int i;

    for (i=0; i<n; i++) v[i]=0.0;
    for (i=0, p=strtok(p, " "); p&&i<n; p=strtok(NULL, " ")) {
        v[i++]=atof(p)*1E-3;
    }
    return i;
}
/* add antenna parameter -----------------------------------------------------*/
static void addpcv(const pcv_t* pcv, pcvs_t* pcvs)
{
    pcv_t* pcvs_pcv;

    if (pcvs->nmax<=pcvs->n) {
        pcvs->nmax+=256;
        if (!(pcvs_pcv=(pcv_t*)realloc(pcvs->pcv, sizeof(pcv_t)*pcvs->nmax))) {
            trace(1, "addpcv: memory allocation error\n");
            free(pcvs->pcv); pcvs->pcv=NULL; pcvs->n=pcvs->nmax=0;
            return;
        }
        pcvs->pcv=pcvs_pcv;
    }
    pcvs->pcv[pcvs->n++]=*pcv;
}
/* read ngs antenna parameter file -------------------------------------------*/
static int readngspcv(const char* file, pcvs_t* pcvs)
{
    FILE* fp;
    static const pcv_t pcv0={ 0 };
    pcv_t pcv;
    double neu[3];
    int n=0;
    char buff[256];

    if (!(fp=fopen(file, "r"))) {
        trace(2, "ngs pcv file open error: %s\n", file);
        return 0;
    }
    while (fgets(buff, sizeof(buff), fp)) {

        if (strlen(buff)>=62&&buff[61]=='|') continue;

        if (buff[0] != ' ') n=0; /* start line */
        if (++n==1) {
            pcv=pcv0;
            strncpy(pcv.type, buff, 61); pcv.type[61]='\0';
        }
        else if (n==2) {
            if (decodef(buff, 3, neu)<3) continue;
            pcv.off[0][0]=neu[1];
            pcv.off[0][1]=neu[0];
            pcv.off[0][2]=neu[2];
        }
        else if (n==3) decodef(buff, 10, pcv.var[0]);
        else if (n==4) decodef(buff, 9, pcv.var[0]+10);
        else if (n==5) {
            if (decodef(buff, 3, neu)<3) continue;;
            pcv.off[1][0]=neu[1];
            pcv.off[1][1]=neu[0];
            pcv.off[1][2]=neu[2];
        }
        else if (n==6) decodef(buff, 10, pcv.var[1]);
        else if (n==7) {
            decodef(buff, 9, pcv.var[1]+10);
            addpcv(&pcv, pcvs);
        }
    }
    fclose(fp);

    return 1;
}
/* read antex file ----------------------------------------------------------*/
static int readantex(const char* file, pcvs_t* pcvs)
{
    FILE* fp;
    static const pcv_t pcv0={ 0 };
    pcv_t pcv;
    double neu[3];
    int i, f, freq=0, state=0, freqs[]={ 1,2,5,6,7,8,0 };
    char buff[256];

    trace(3, "readantex: file=%s\n", file);

    if (!(fp=fopen(file, "r"))) {
        trace(2, "antex pcv file open error: %s\n", file);
        return 0;
    }
    while (fgets(buff, sizeof(buff), fp)) {

        if (strlen(buff)<60||strstr(buff+60, "COMMENT")) continue;

        if (strstr(buff+60, "START OF ANTENNA")) {
            pcv=pcv0;
            state=1;
        }
        if (strstr(buff+60, "END OF ANTENNA")) {
            addpcv(&pcv, pcvs);
            state=0;
        }
        if (!state) continue;

        if (strstr(buff+60, "TYPE / SERIAL NO")) {
            strncpy(pcv.type, buff, 20); pcv.type[20]='\0';
            strncpy(pcv.code, buff+20, 20); pcv.code[20]='\0';
            if (!strncmp(pcv.code+3, "        ", 8)) {
                pcv.sat=satid2no(pcv.code);
            }
        }
        else if (strstr(buff+60, "VALID FROM")) {
            if (!str2time(buff, 0, 43, &pcv.ts)) continue;
        }
        else if (strstr(buff+60, "VALID UNTIL")) {
            if (!str2time(buff, 0, 43, &pcv.te)) continue;
        }
        else if (strstr(buff+60, "START OF FREQUENCY")) {
            if (sscanf(buff+4, "%d", &f)<1) continue;
            for (i=0; i<NFREQ; i++) if (freqs[i]==f) break;
            if (i<NFREQ) freq=i+1;
        }
        else if (strstr(buff+60, "END OF FREQUENCY")) {
            freq=0;
        }
        else if (strstr(buff+60, "NORTH / EAST / UP")) {
            if (freq<1||NFREQ<freq) continue;
            if (decodef(buff, 3, neu)<3) continue;
            pcv.off[freq-1][0]=neu[pcv.sat?0:1]; /* x or e */
            pcv.off[freq-1][1]=neu[pcv.sat?1:0]; /* y or n */
            pcv.off[freq-1][2]=neu[2];           /* z or u */
        }
        else if (strstr(buff, "NOAZI")) {
            if (freq<1||NFREQ<freq) continue;
            if ((i=decodef(buff+8, 19, pcv.var[freq-1]))<=0) continue;
            for (; i<19; i++) pcv.var[freq-1][i]=pcv.var[freq-1][i-1];
        }
    }
    fclose(fp);

    return 1;
}
/* read antenna parameters ------------------------------------------------------
* read antenna parameters
* args  :char   *file       I   antenna parameter file (antex)
*          pcvs_t *pcvs       IO  antenna parameters
* return:status (1:ok,0:file open error)
* notes :file with the externsion .atx or .ATX is recognized as antex
*          file except for antex is recognized ngs antenna parameters
*          see reference [3]
*          only support non-azimuth-depedent parameters
*-----------------------------------------------------------------------------*/
extern int readpcv(const char* file, pcvs_t* pcvs)
{
    pcv_t* pcv;
    char* ext;
    int i, stat;

    trace(3, "readpcv: file=%s\n", file);

    if (!(ext=strrchr(file, '.'))) ext="";

    if (!strcmp(ext, ".atx")||!strcmp(ext, ".ATX")) {
        stat=readantex(file, pcvs);
    }
    else {
        stat=readngspcv(file, pcvs);
    }
    for (i=0; i<pcvs->n; i++) {
        pcv=pcvs->pcv+i;
        trace(4, "sat=%2d type=%20s code=%s off=%8.4f %8.4f %8.4f  %8.4f %8.4f %8.4f\n",
            pcv->sat, pcv->type, pcv->code, pcv->off[0][0], pcv->off[0][1],
            pcv->off[0][2], pcv->off[1][0], pcv->off[1][1], pcv->off[1][2]);
    }
    return stat;
}