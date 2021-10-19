
#include "io.h"

/* read DCB parameters file --------------------------------------------------*/
static int readdcbf(const char *file, nav_t *nav, const sta_t *sta)
{
    FILE *fp;
    double cbias;
    char buff[256],str1[32],str2[32]="";
    int i,j,sat,type=0;

    trace(3,"readdcbf: file=%s\n",file);

    if (!(fp=fopen(file,"r"))) {
        trace(2,"dcb parameters file open error: %s\n",file);
        return 0;
    }
    while (fgets(buff,sizeof(buff),fp)) {

        if      (strstr(buff,"DIFFERENTIAL (P1-P2) CODE BIASES")) type=1;
        else if (strstr(buff,"DIFFERENTIAL (P1-C1) CODE BIASES")) type=2;
        else if (strstr(buff,"DIFFERENTIAL (P2-C2) CODE BIASES")) type=3;

        if (!type||sscanf(buff,"%s %s",str1,str2)<1) continue;

        if ((cbias=str2num(buff,26,9))==0.0) continue;

        if (sta&&(!strcmp(str1,"G")||!strcmp(str1,"R"))) { /* receiver DCB */
            for (i=0;i<MAXRCV;i++) {
                if (!strcmp(sta[i].name,str2)) break;
            }
            if (i<MAXRCV) {
                j=!strcmp(str1,"G")?0:1;
                nav->rbias[i][j][type-1]=cbias*1E-9*CLIGHT; /* ns -> m */
            }
        }
        else if ((sat=satid2no(str1))) { /* satellite dcb */
            nav->cbias[sat-1][type-1]=cbias*1E-9*CLIGHT; /* ns -> m */
        }
    }
    fclose(fp);

    return 1;
}
/* read DCB parameters ---------------------------------------------------------
* read differential code bias (DCB) parameters
* args   : char   *file       I   DCB parameters file (wild-card * expanded)
*          nav_t  *nav        IO  navigation data
*          sta_t  *sta        I   station info data to inport receiver DCB
*                                 (NULL: no use)
* return : status (1:ok,0:error)
* notes  : currently only support P1-P2, P1-C1, P2-C2, bias in DCB file
*-----------------------------------------------------------------------------*/
extern int readdcb(const char *file, nav_t *nav, const sta_t *sta)
{
    int i,j,n;
    char *efiles[MAXEXFILE]={0};

    trace(3,"readdcb : file=%s\n",file);

    for (i=0;i<MAXSAT;i++) for (j=0;j<3;j++) {
        nav->cbias[i][j]=0.0;
    }
    for (i=0;i<MAXEXFILE;i++) {
        if (!(efiles[i]=(char *)malloc(1024))) {
            for (i--;i>=0;i--) free(efiles[i]);
            return 0;
        }
    }
    n=expath(file,efiles,MAXEXFILE);

    for (i=0;i<n;i++) {
        readdcbf(efiles[i],nav,sta);
    }
    for (i=0;i<MAXEXFILE;i++) free(efiles[i]);

    return 1;
}