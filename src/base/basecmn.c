
#include "base.h"

/* constants -----------------------------------------------------------------*/
const double chisqr[100]={      /* chi-sqr(n) (alpha=0.001) */
    10.8,13.8,16.3,18.5,20.5,22.5,24.3,26.1,27.9,29.6,
    31.3,32.9,34.5,36.1,37.7,39.3,40.8,42.3,43.8,45.3,
    46.8,48.3,49.7,51.2,52.6,54.1,55.5,56.9,58.3,59.7,
    61.1,62.5,63.9,65.2,66.6,68.0,69.3,70.7,72.1,73.4,
    74.7,76.0,77.3,78.6,80.0,81.3,82.6,84.0,85.4,86.7,
    88.0,89.3,90.6,91.9,93.3,94.7,96.0,97.4,98.7,100 ,
    101 ,102 ,103 ,104 ,105 ,107 ,108 ,109 ,110 ,112 ,
    113 ,114 ,115 ,116 ,118 ,119 ,120 ,122 ,123 ,125 ,
    126 ,127 ,128 ,129 ,131 ,132 ,133 ,134 ,135 ,137 ,
    138 ,139 ,140 ,142 ,143 ,144 ,145 ,147 ,148 ,149
};
const prcopt_t prcopt_default={ /* defaults processing options */
    {0},{0},0.0,0.0,            /* ts,te,ti,tu */
    {0},1,                      /* network processing time */
    PMODE_SINGLE,0,2,SYS_GPS,   /* mode,soltype,nf,navsys */
    15.0*D2R,{{0,0}},           /* elmin,snrmask */
    EPHOPT_PREC,ARMODE_CONT,    /* sateph,modear */
    IONOOPT_EST,TROPOPT_EST,    /* estion,esttrop */
    0.0,                        /* elmaskar */
    30.0,                       /* maxgdop */
    0,                          /* rovpos */
    {100.0,100.0},              /* eratio[] */
    {100.0,0.003,0.003,0.0,1.0},/* err[] */
    {0.0,1E-7,100.0,0.0,0.0,0.0,100.0,100.0},     /* prn */
    {3.0,0.9999,0.25,0.1,0.05}, /* thresar */
    {0},                        /* ru */
    {"",""},{{0}},              /* anttype, antdel */
    {{0}},{{0}},                /* pcvr,exsats */
    {"",""},{{0}},{""},         /* rnxopt, odisp, ppopt */
    1,0,1,0,                    /* ircv, inh, estclk, datumType */
    0,0,0                       /* rcb, scb, dck */
};

const solopt_t solopt_default={ /* defaults solution output options */
    SOLF_XYZ,TIMES_GPST,1,2,    /* posf,times,timef,timeu */
    0,1,1,0,0,0,0,              /* degf,outhead,outopt,outvel,datum,height,geoid */
    0,0,0,                      /* solstatic,sstat,trace */
    {0.0,0.0},                  /* nmeaintv */
    " ","",0.0                  /* separator/program name */
};

/* execute command -------------------------------------------------------------
* execute command line by operating system shell
* args   : char   *cmd      I   command line
* return : execution status (0:ok,0>:error)
*-----------------------------------------------------------------------------*/
extern int execcmd(const char *cmd)
{
#ifdef WIN32
    PROCESS_INFORMATION info;
    STARTUPINFO si={0};
    DWORD stat;
    char cmds[1024];

    trace(3,"execcmd: cmd=%s\n",cmd);

    si.cb=sizeof(si);
    sprintf(cmds,"cmd /c %s",cmd);
    if (!CreateProcess(NULL,(LPTSTR)cmds,NULL,NULL,FALSE,CREATE_NO_WINDOW,NULL,
        NULL,&si,&info)) return -1;
    WaitForSingleObject(info.hProcess,INFINITE);
    if (!GetExitCodeProcess(info.hProcess,&stat)) stat=-1;
    CloseHandle(info.hProcess);
    CloseHandle(info.hThread);
    return (int)stat;
#else
    trace(3,"execcmd: cmd=%s\n",cmd);

    return system(cmd);
#endif
}
/* uncompress file -------------------------------------------------------------
* uncompress (uncompress/unzip/uncompact hatanaka-compression/tar) file
* args   : char   *file     I   input file
*          char   *uncfile  O   uncompressed file
* return : status (-1:error,0:not compressed file,1:uncompress completed)
* note   : creates uncompressed file in tempolary directory
*          gzip, tar and crx2rnx commands have to be installed in commands path
*-----------------------------------------------------------------------------*/
extern int rtk_uncompress(const char *file, char *uncfile)
{
    int stat=0;
    char *p,cmd[64+2048]="",tmpfile[1024]="",buff[1024],*fname,*dir="";

    trace(3,"rtk_uncompress: file=%s\n",file);

    strcpy(tmpfile,file);
    if (!(p=strrchr(tmpfile,'.'))) return 0;

    /* uncompress by gzip */
    if (!strcmp(p,".z"  )||!strcmp(p,".Z"  )||
        !strcmp(p,".gz" )||!strcmp(p,".GZ" )||
        !strcmp(p,".zip")||!strcmp(p,".ZIP")) {

        strcpy(uncfile,tmpfile); uncfile[p-tmpfile]='\0';
        sprintf(cmd,"gzip -f -d -c \"%s\" > \"%s\"",tmpfile,uncfile);

        if (execcmd(cmd)) {
            remove(uncfile);
            return -1;
        }
        strcpy(tmpfile,uncfile);
        stat=1;
    }
    /* extract tar file */
    if ((p=strrchr(tmpfile,'.'))&&!strcmp(p,".tar")) {

        strcpy(uncfile,tmpfile); uncfile[p-tmpfile]='\0';
        strcpy(buff,tmpfile);
        fname=buff;
#ifdef WIN32
        if ((p=strrchr(buff,'\\'))) {
            *p='\0'; dir=fname; fname=p+1;
        }
        sprintf(cmd,"set PATH=%%CD%%;%%PATH%% & cd /D \"%s\" & tar -xf \"%s\"",
            dir,fname);
#else
        if ((p=strrchr(buff,'/'))) {
            *p='\0'; dir=fname; fname=p+1;
        }
        sprintf(cmd,"tar -C \"%s\" -xf \"%s\"",dir,tmpfile);
#endif
        if (execcmd(cmd)) {
            if (stat) remove(tmpfile);
            return -1;
        }
        if (stat) remove(tmpfile);
        stat=1;
    }
    /* extract hatanaka-compressed file by cnx2rnx */
    else if ((p=strrchr(tmpfile,'.'))&&
        ((strlen(p)>3&&(*(p+3)=='d'||*(p+3)=='D'))||
            !strcmp(p,".crx")||!strcmp(p,".CRX"))) {

        strcpy(uncfile,tmpfile);
        uncfile[p-tmpfile+3]=*(p+3)=='D'?'O':'o';
        sprintf(cmd,"crx2rnx < \"%s\" > \"%s\"",tmpfile,uncfile);

        if (execcmd(cmd)) {
            remove(uncfile);
            if (stat) remove(tmpfile);
            return -1;
        }
        if (stat) remove(tmpfile);
        stat=1;
    }
    trace(3,"rtk_uncompress: stat=%d\n",stat);
    return stat;
}