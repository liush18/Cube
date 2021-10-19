
#include "net.h"

static FILE *fpDatum=NULL;    /* file pointer of ambiguity datum */

/* open ambiguity datum file -------------------------------------------------*/
extern int openDtmLog(const char* file)
{
    createdir(file);    /* by ls */
    if (!(fpDatum=fopen(file, "w"))) {
        trace(2,"datum_open: file open error path=%s\n",file);
        return 0;
    }

    return 1;
}
/* close ambiguity datum file ------------------------------------------------*/
extern void closeDtmLog(void)
{
    if (fpDatum) fclose(fpDatum);
    fpDatum=NULL;
}
extern void logDatumRcvPos(rcv_t *rcv, const int n)
{
    double pos[3];
    char cst_str[32]="RCV/POSITION";
    int i;

    fprintf(fpDatum,"+%s\n",cst_str);
    for (i=0;i<n&&i<MAXRCV;i++) {
        ecef2pos(rcv[i].sta.pos,pos);
        fprintf(fpDatum,"%-4s %8.3f %8.3f %14.3f\n",
            rcv[i].sta.name,pos[0]*R2D,pos[1]*R2D,pos[2]);
    }
    fprintf(fpDatum,"-%s\n",cst_str);
    fflush(fpDatum);
}
/* output ambiguity datum imformation ----------------------------------------*/
extern void logDtm(gtime_t time, rcv_t *rcv, const int n)
{
    double pos[3],*rs;
    char str[32],id[4],cst_str1[32]="SAT/POSITION",cst_str2[32]="DATUM";
    int i,j,sat;
    
    if (!fpDatum) return;

    time2str(time,str,2);
    fprintf(fpDatum,"> %s\n",str);

    /* receiver and satellite position */
    fprintf(fpDatum,"+%s\n",cst_str1);
    rs=zeros(3,MAXSAT);
    for (i=0;i<n&&i<MAXRCV;i++) {
        if (rcv[i].nobs==0) continue;
        for (j=0;j<rcv[i].nobs&&j<MAXOBS;j++) {
            sat=rcv[i].obs[j].sat;
            if (rcv[i].datum[sat-1]==OBS_REJ) continue;
            if (norm(rs+(sat-1)*3,3)<1.0&&norm(rcv[i].rs+j*3,3)>1.0)
                matcpy(rs+(sat-1)*3,rcv[i].rs+j*6,3,1);
        }
    }
    for (i=0;i<MAXSAT;i++) {
        if (norm(rs+i*3,3)<1.0) continue;
        satno2id(i+1,id);
        ecef2pos(rs+i*3,pos);
        fprintf(fpDatum,"%-4s %8.3f %8.3f %14.3f\n",
            id,pos[0]*R2D,pos[1]*R2D,pos[2]);
    }
    fprintf(fpDatum,"-%s\n",cst_str1);

    /* ambiguity datum (0:not observed, 1:objected, 2:observed, 3:datum) */
    fprintf(fpDatum,"+%s\n",cst_str2);
    fprintf(fpDatum,"%4s","");
    for (i=0;i<MAXSAT;i++) {
        satno2id(i+1,id);
        fprintf(fpDatum,"%4s",id);
    }
    fprintf(fpDatum,"\n");
    for (i=0;i<n&&i<MAXRCV;i++) {
        if (rcv[i].nobs==0) continue;
        fprintf(fpDatum,"%4s",rcv[i].sta.name);
        for (j=0;j<MAXSAT;j++) fprintf(fpDatum,"%4d",rcv[i].datum[j]);
        fprintf(fpDatum,"\n");
    }
    fprintf(fpDatum,"-%s\n",cst_str2);
    fflush(fpDatum);
}