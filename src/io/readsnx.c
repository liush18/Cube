
#if 0

#include "io.h"


/* get coordinate for stations -----------------------------------------------*/
extern int readsnx(const char *file, rcv_t *rcv, const int n)
{
	FILE *fp;
	char buff[1024],stalist[MAXRCV*5+1],code[5],*p;
	int i,flag=0;

	trace(3,"readsnx:file=%s\n",file);

	if (!(fp=fopen(file,"r"))) {
		trace(2,"snx file open error: %s\n",file);
		return 0;
	}

	p=stalist;
	for (i=0;i<n;i++) {
		rcv[i].rr[0]=rcv[i].rr[1]=rcv[i].rr[2]=0.0;
		strncpy(p,rcv[i].name,4);
		p+=4;
		sprintf(p++,"%c",' ');
	}
		
	while (fgets(buff,sizeof(buff),fp)) {
		if (strstr(buff,"+SOLUTION/ESTIMATE")) {
			flag=1;
			continue;
		}
		if (strstr(buff,"-SOLUTION/ESTIMATE")) break;

		strncpy(code,&buff[14],4);code[4]='\0';
		if (flag&&strstr(stalist,code)) {
			for (i=0;i<n;i++) {
				if (rcv[i].excrcv) continue;
				if (!strncmp(code,rcv[i].name,4)&&!strncmp(&buff[7],"STAX",4)) {
					rcv[i].rr[0]=str2num(buff,47,21);
					break;
				}	
				else if (!strncmp(code,rcv[i].name,4)&&!strncmp(&buff[7],"STAY",4)) {
					rcv[i].rr[1]=str2num(buff,47,21);
					break;
				}	
				else if (!strncmp(code,rcv[i].name,4)&&!strncmp(&buff[7],"STAZ",4)) {
					rcv[i].rr[2]=str2num(buff,47,21);
					break;
				}	
			}

		}
	}

	for (i=0;i<n;i++) {
		if (rcv[i].excrcv) continue;
		if (rcv[i].rr[0]*rcv[i].rr[1]*rcv[i].rr[2]==0.0) {
			rcv[i].excrcv=1;
			trace(2,"no snx coordinates: sta=%s\n",rcv[i].name);
		}	
	}

	fclose(fp);
	return 1;
}

#if 0

/* get coordinate for stations -----------------------------------------------*/
extern int readsnx(const char *file,sta_t *sta,const int n)
{
	FILE *fp;
	char buff[1024],stalist[MAXRCV*5+1],code[5],*p;
	int i,flag=0;

	trace(3,"readsnx:file=%s\n",file);

	if (!(fp=fopen(file,"r"))) {
		trace(2,"snx file open error: %s\n",file);
		return 0;
	}

	p=stalist;
	for (i=0;i<n;i++) {
		sta[i].pos[0]=sta[i].pos[1]=sta[i].pos[2]=0.0;
		strncpy(p,sta[i].name,4);
		p+=4;
		sprintf(p++,"%c",' ');
	}

	while (fgets(buff,sizeof(buff),fp)) {
		if (strstr(buff,"+SOLUTION/ESTIMATE")) {
			flag=1;
			continue;
		}
		if (strstr(buff,"-SOLUTION/ESTIMATE")) break;

		strncpy(code,&buff[14],4);code[4]='\0';
		if (flag&&strstr(stalist,code)) {
			for (i=0;i<n;i++) {
				if (sta[i].exc) continue;
				if (!strncmp(code,sta[i].name,4)&&!strncmp(&buff[7],"STAX",4)) {
					sta[i].pos[0]=str2num(buff,47,21);
					break;
				}	
				else if (!strncmp(code,sta[i].name,4)&&!strncmp(&buff[7],"STAY",4)) {
					sta[i].pos[1]=str2num(buff,47,21);
					break;
				}	
				else if (!strncmp(code,sta[i].name,4)&&!strncmp(&buff[7],"STAZ",4)) {
					sta[i].pos[2]=str2num(buff,47,21);
					break;
				}	
			}

		}
	}

	for (i=0;i<n;i++) {
		if (sta[i].exc) continue;
		if (sta[i].pos[0]*sta[i].pos[1]*sta[i].pos[2]==0.0) {
			sta[i].exc=1;
			trace(2,"no snx coordinates: sta=%s\n",sta[i].name);
		}	
	}

	fclose(fp);
	return 1;
}
#endif


#endif