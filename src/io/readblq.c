
#include "io.h"

/* read blq record -----------------------------------------------------------*/
static int readblqrecord(FILE* fp, double* odisp)
{
    double v[11];
    char buff[256];
    int i, n=0;

    while (fgets(buff, sizeof(buff), fp)) {
        if (!strncmp(buff, "$$", 2)) continue;
        if (sscanf(buff, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            v, v+1, v+2, v+3, v+4, v+5, v+6, v+7, v+8, v+9, v+10)<11) continue;
        for (i=0; i<11; i++) odisp[n+i*6]=v[i];
        if (++n==6) return 1;
    }
    return 0;
}
/* read blq ocean tide loading parameters --------------------------------------
* read blq ocean tide loading parameters
* args  :char   *file       I   BLQ ocean tide loading parameter file
*          char   *sta        I   station name
*          double *odisp      O   ocean tide loading parameters
* return:status (1:ok,0:file open error)
*-----------------------------------------------------------------------------*/
extern int readblq(const char* file, const char* sta, double* odisp)
{
    FILE* fp;
    char buff[256], staname[32]="", name[32],*p;

    /* station name to upper case */
    sscanf(sta, "%16s", staname);
    for (p=staname; (*p=(char)toupper((int)(*p))); p++);

    if (!(fp=fopen(file, "r"))) {
        trace(2, "blq file open error: file=%s\n", file);
        return 0;
    }
    while (fgets(buff, sizeof(buff), fp)) {
        if (!strncmp(buff, "$$", 2)||strlen(buff)<2) continue;

        if (sscanf(buff+2, "%16s", name)<1) continue;
        for (p=name; (*p=(char)toupper((int)(*p))); p++);
        if (strcmp(name, staname)) continue;

        /* read blq record */
        if (readblqrecord(fp, odisp)) {
            fclose(fp);
            return 1;
        }
    }
    fclose(fp);
    trace(2, "no otl parameters: sta=%s file=%s\n",sta,file);
    return 0;
}
/* read ocean tide loading parameters ----------------------------------------*/
extern void readotl_pos(prcopt_t *popt, const char *file, const sta_t *sta)
{
    int i,mode=PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_FIXED;

    for (i=0;i<(mode?2:1);i++) {
        readblq(file,sta[i].name,popt->odisp[i]);
    }
}