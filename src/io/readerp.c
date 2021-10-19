
#include "io.h"


/* read earth rotation parameters ----------------------------------------------
* read earth rotation parameters
* args  :char   *file       I   IGS ERP file (IGS ERP ver.2)
*          erp_t  *erp        O   earth rotation parameters
* return:status (1:ok,0:file open error)
*-----------------------------------------------------------------------------*/
extern int readerp(const char* file, erp_t* erp)
{
    FILE* fp;
    erpd_t* erp_data;
    double v[14]={ 0 };
    char buff[256];

    trace(3, "readerp: file=%s\n", file);

    if (!(fp=fopen(file, "r"))) {
        trace(2, "erp file open error: file=%s\n", file);
        return 0;
    }
    while (fgets(buff, sizeof(buff), fp)) {
        if (sscanf(buff, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            v, v+1, v+2, v+3, v+4, v+5, v+6, v+7, v+8, v+9, v+10, v+11, v+12, v+13)<5) {
            continue;
        }
        if (erp->n>=erp->nmax) {
            erp->nmax=erp->nmax<=0?128:erp->nmax*2;
            erp_data=(erpd_t*)realloc(erp->data, sizeof(erpd_t)*erp->nmax);
            if (!erp_data) {
                free(erp->data); erp->data=NULL; erp->n=erp->nmax=0;
                fclose(fp);
                return 0;
            }
            erp->data=erp_data;
        }
        erp->data[erp->n].mjd=v[0];
        erp->data[erp->n].xp=v[1]*1E-6*AS2R;
        erp->data[erp->n].yp=v[2]*1E-6*AS2R;
        erp->data[erp->n].ut1_utc=v[3]*1E-7;
        erp->data[erp->n].lod=v[4]*1E-7;
        erp->data[erp->n].xpr=v[12]*1E-6*AS2R;
        erp->data[erp->n++].ypr=v[13]*1E-6*AS2R;
    }
    fclose(fp);
    return 1;
}