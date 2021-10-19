
#include "base.h"

/* string to number ------------------------------------------------------------
* convert substring in string to number
* args  :  char   *s        I   string ("... nnn.nnn ...")
*          int    i,n       I   substring position and width
* return:converted number (0.0:error)
*-----------------------------------------------------------------------------*/
extern double str2num(const char *s, int i, int n)
{
	double value;
	char str[256], *p=str;

	if (i<0||(int)strlen(s)<i||(int)sizeof(str)-1<n) return 0.0;
	for (s+=i; *s&&--n>=0; s++) *p++=*s=='d'||*s=='D'?'E':*s;
	*p='\0';
	return sscanf(str, "%lf", &value)==1?value:0.0;
}
/* expand file path ------------------------------------------------------------
* expand file path with wild-card (*) in file
* args  :char   *path     I   file path to expand (captal insensitive)
*          char   *paths    O   expanded file paths
*          int    nmax      I   max number of expanded file paths
* return:number of expanded file paths
* notes :the order of expanded files is alphabetical order
*-----------------------------------------------------------------------------*/
extern int expath(const char *path, char *paths[], int nmax)
{
	int i, j, n=0;
	char tmp[1024];
#ifdef WIN32
	WIN32_FIND_DATA file;
	HANDLE h;
	char dir[1024]="", *p;

	trace(3, "expath :path=%s nmax=%d\n", path, nmax);

	if ((p=strrchr(path, '\\'))) {
		strncpy(dir, path, p-path+1); dir[p-path+1]='\0';
	}
	if ((h=FindFirstFile((LPCTSTR)path, &file))==INVALID_HANDLE_VALUE) {
		strcpy(paths[0], path);
		return 1;
	}
	sprintf(paths[n++], "%s%s", dir, file.cFileName);
	while (FindNextFile(h, &file)&&n<nmax) {
		if (file.dwFileAttributes&FILE_ATTRIBUTE_DIRECTORY) continue;
		sprintf(paths[n++], "%s%s", dir, file.cFileName);
	}
	FindClose(h);
#else
	struct dirent *d;
	DIR *dp;
	const char *file=path;
	char dir[1024]="", s1[1024], s2[1024], *p, *q, *r;

	trace(3, "expath :path=%s nmax=%d\n", path, nmax);

	if ((p=strrchr(path, '/'))||(p=strrchr(path, '\\'))) {
		file=p+1; strncpy(dir, path, p-path+1); dir[p-path+1]='\0';
	}
	if (!(dp=opendir(*dir?dir:"."))) return 0;
	while ((d=readdir(dp))) {
		if (*(d->d_name)=='.') continue;
		sprintf(s1, "^%s$", d->d_name);
		sprintf(s2, "^%s$", file);
		for (p=s1; *p; p++) *p=(char)tolower((int)*p);
		for (p=s2; *p; p++) *p=(char)tolower((int)*p);

		for (p=s1, q=strtok_r(s2, "*", &r); q; q=strtok_r(NULL, "*", &r)) {
			if ((p=strstr(p, q))) p+=strlen(q); else break;
		}
		if (p&&n<nmax) sprintf(paths[n++], "%s%s", dir, d->d_name);
	}
	closedir(dp);
#endif
	/* sort paths in alphabetical order */
	for (i=0; i<n-1; i++) {
		for (j=i+1; j<n; j++) {
			if (strcmp(paths[i], paths[j])>0) {
				strcpy(tmp, paths[i]);
				strcpy(paths[i], paths[j]);
				strcpy(paths[j], tmp);
			}
		}
	}
	for (i=0; i<n; i++) trace(3, "expath :file=%s\n", paths[i]);

	return n;
}
/* create directory ------------------------------------------------------------
* create directory if not exist
* args  :char   *path     I   file path to be saved
* return:none
* notes :not recursive. only one level
*-----------------------------------------------------------------------------*/
extern void createdir(const char *path)
{
	char buff[1024], *p;

	tracet(3, "createdir: path=%s\n", path);

	strcpy(buff, path);
	if (!(p=strrchr(buff, FILEPATHSEP))) return;
	*p='\0';

#ifdef WIN32
	CreateDirectory(buff, NULL);
#else
	mkdir(buff, 0777);
#endif
}
/* set string without tail space ---------------------------------------------*/
extern void setstr(char *dst, const char *src, int n)
{
	char *p=dst;
	const char *q=src;
	while (*q&&q<src+n) *p++=*q++;
	*p--='\0';
	while (p>=dst&&*p==' ') *p--='\0';
}
/* convert string to upper ---------------------------------------------------*/
extern void strup(char *p)
{
	for (; (*p=(char)toupper((int)(*p))); p++);
}
/* discard space characters at tail ------------------------------------------*/
extern void chop(char *str)
{
	char *p;
	if ((p=strchr(str,'#'))) *p='\0'; /* comment */
	for (p=str+strlen(str)-1;p>=str&&!isgraph((int)*p);p--) *p='\0';
}
/* get o file path from station list  ----------------------------------------*/
extern int getOFilePath(const char *file, char *files[])
{
	FILE *fp;
	char buff[1024];
	int i=0;

	if (!(fp=fopen(file,"r"))) {
		trace(1,"file open error: %s\n",file);
		return 0;
	}
	while (fgets(buff,1024,fp)) {
		if (strlen(buff)<5) continue;
		buff[strlen(buff)-1]='\0';
		sprintf(files[i],buff);
		i++;
	}
	fclose(fp);
	return i;
}
/* show message --------------------------------------------------------------*/
extern int showmsg(char* format, ...)
{
	va_list arg;
	va_start(arg, format); vfprintf(stderr, format, arg); va_end(arg);
	fprintf(stderr, "\r");
	//fprintf(stderr, "\n");
	return 0;
}
/* replace string ------------------------------------------------------------*/
static int repstr(char *str, const char *pat, const char *rep)
{
	int len=(int)strlen(pat);
	char buff[1024],*p,*q,*r;

	for (p=str,r=buff;*p;p=q+len) {
		if (!(q=strstr(p,pat))) break;
		strncpy(r,p,q-p);
		r+=q-p;
		r+=sprintf(r,"%s",rep);
	}
	if (p<=str) return 0;
	strcpy(r,p);
	strcpy(str,buff);
	return 1;
}
/* replace keywords in file path -----------------------------------------------
* replace keywords in file path with date, time, rover and base station id
* args   : char   *path     I   file path (see below)
*          char   *rpath    O   file path in which keywords replaced (see below)
*          gtime_t time     I   time (gpst)  (time.time==0: not replaced)
*          char   *rov      I   rover id string        ("": not replaced)
*          char   *base     I   base station id string ("": not replaced)
* return : status (1:keywords replaced, 0:no valid keyword in the path,
*                  -1:no valid time)
* notes  : the following keywords in path are replaced by date, time and name
*              %Y -> yyyy : year (4 digits) (1900-2099)
*              %y -> yy   : year (2 digits) (00-99)
*              %m -> mm   : month           (01-12)
*              %d -> dd   : day of month    (01-31)
*              %h -> hh   : hours           (00-23)
*              %M -> mm   : minutes         (00-59)
*              %S -> ss   : seconds         (00-59)
*              %n -> ddd  : day of year     (001-366)
*              %W -> wwww : gps week        (0001-9999)
*              %D -> d    : day of gps week (0-6)
*              %H -> h    : hour code       (a=0,b=1,c=2,...,x=23)
*              %ha-> hh   : 3 hours         (00,03,06,...,21)
*              %hb-> hh   : 6 hours         (00,06,12,18)
*              %hc-> hh   : 12 hours        (00,12)
*              %t -> mm   : 15 minutes      (00,15,30,45)
*-----------------------------------------------------------------------------*/
extern int reppath(const char *path, char *rpath, gtime_t time)
{
	double ep[6],ep0[6]={2000,1,1,0,0,0};
	int week,dow,doy,stat=0;
	char rep[64];

	strcpy(rpath,path);

	if (!strstr(rpath,"%")) return 0;
	if (time.time!=0) {
		time2epoch(time,ep);
		ep0[0]=ep[0]; /* replace year */
		dow=(int)floor(time2gpst(time,&week)/86400.0); /* day of gps week */
		doy=(int)floor(timediff(time,epoch2time(ep0))/86400.0)+1; /* day of year */
		sprintf(rep,"%02d",  ((int)ep[3]/3)*3);   stat|=repstr(rpath,"%ha",rep);
		sprintf(rep,"%02d",  ((int)ep[3]/6)*6);   stat|=repstr(rpath,"%hb",rep);
		sprintf(rep,"%02d",  ((int)ep[3]/12)*12); stat|=repstr(rpath,"%hc",rep);
		sprintf(rep,"%04.0f",ep[0]);              stat|=repstr(rpath,"%Y",rep);
		sprintf(rep,"%02.0f",fmod(ep[0],100.0));  stat|=repstr(rpath,"%y",rep);
		sprintf(rep,"%02.0f",ep[1]);              stat|=repstr(rpath,"%m",rep);
		sprintf(rep,"%02.0f",ep[2]);              stat|=repstr(rpath,"%d",rep);
		sprintf(rep,"%02.0f",ep[3]);              stat|=repstr(rpath,"%h",rep);
		sprintf(rep,"%02.0f",ep[4]);              stat|=repstr(rpath,"%M",rep);
		sprintf(rep,"%02.0f",floor(ep[5]));       stat|=repstr(rpath,"%S",rep);
		sprintf(rep,"%03d",  doy);                stat|=repstr(rpath,"%n",rep);
		sprintf(rep,"%04d",  week);               stat|=repstr(rpath,"%W",rep);
		sprintf(rep,"%d",    dow);                stat|=repstr(rpath,"%D",rep);
		sprintf(rep,"%c",    'a'+(int)ep[3]);     stat|=repstr(rpath,"%H",rep);
		sprintf(rep,"%02d",  ((int)ep[4]/15)*15); stat|=repstr(rpath,"%t",rep);
	}
	else if (strstr(rpath,"%ha")||strstr(rpath,"%hb")||strstr(rpath,"%hc")||
		strstr(rpath,"%Y" )||strstr(rpath,"%y" )||strstr(rpath,"%m" )||
		strstr(rpath,"%d" )||strstr(rpath,"%h" )||strstr(rpath,"%M" )||
		strstr(rpath,"%S" )||strstr(rpath,"%n" )||strstr(rpath,"%W" )||
		strstr(rpath,"%D" )||strstr(rpath,"%H" )||strstr(rpath,"%t" )) {
		return -1; /* no valid time */
	}
	return stat;
}