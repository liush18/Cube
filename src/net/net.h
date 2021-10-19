
#ifndef _NET_H
#define _NET_H

#include "corr/corr.h"
#include "io/io.h"
#include "qc/qc.h"
#include "sat/sat.h"


#define OBS_NONE    0       /* no obs data */
#define OBS_REJ     1       /* reject obs data */
#define OBS_USE     2       /* obs data in use */
#define OBS_SSC     3       /* ambiguity datum for separating satellites' clocks */
#define OBS_SRC     4       /* ambiguity datum for separating receivers' clocks */
#define OBS_FIX		5		/* fixed ambiguity */

#define MAXSATOUT	1		/* max outage counter of satellite */
#define MAXRCVOUT	1		/* max outage counter of receiver */

/* constraint type */
#define TYPE_TIME	0		/* time */
#define TYPE_TRP	1		/* troposphere */
#define TYPE_ION	2		/* ionosphere */
#define TYPE_AMB	3		/* ambiguity */
#define	TYPE_SCB	4		/* satellite code bias */
#define TYPE_RCB	5		/* receiver code bias */
#define TYPE_SPB	6		/* satellite phase bias */
#define TYPE_RPB	7		/* receiver phase bias */

/* constraint string flag */
#define CST_STR_TIME	"#TIME"
#define CST_STR_TRP		"#TRP"
#define CST_STR_ION		"#ION"
#define CST_STR_AMB		"#AMB"
#define CST_STR_SCB		"#SCB"
#define CST_STR_RCB		"#RCB"
#define CST_STR_SPB		"#SPB"
#define CST_STR_RPB		"#RPB"

#define DTM_ZMCSAT	0		/* satellite datum */
#define DTM_REFRCV	1		/* receiver datum */

#define FST_EPOCH	0		/* first epoch */
#define MID_EPOCH	1		/* middle epoch */
#define LST_EPOCH	2		/* last epoch */

#define MAX_REJTIMES	3	/* max reject times in filter */
#define CONSFIXED		5	/* continuous fixed times of ambiguity hold */
#define NET_TIMEZONE	8	/* time zone */

typedef struct {		/* constraint data struct */
	int type;			/* 1:sbias,2:trop,3:amb,4:rbias */
	char rcv[8];		/* receiver name with four characters */
	int sat;			/* satellite number */
	int frq;			/* frequency number */
	int trp;			/* trop type, from 0 */
	double value;		/* estimated values for constraint */
}csd_t;

typedef struct {		/* ambiguity datum struct */
	char rcv[8];		/* receiver name with four characters */
	int sat;			/* satellite number */
	int frq;			/* frequency number */
	int datum;			/* datum for ambiguity */
	double value;		/* estimated values for constraint */
}dtm_t;

typedef struct {		/* observation RINEX record */
	double ver;			/* version */
	int sys,tsys;		/* system and time system */
	char tobs[NUMSYS][MAXOBSTYPE][4];	/* type of observation */
	int flag;			/* observation data flag */
	gtime_t time;		/* current obs time in rcv_t obsd_t */
	int n;				/* current obs number in rcv_t obsd_t */
}rnxo_t;

typedef struct {
	gtime_t time;		/* time (GPST) */
	int stat;
	float ratio;			/* AR ratio factor for valiation */
} netsol_t;

typedef struct {
	int nx;				/* number of float states */
	double *x,*P;		/* float states and their covariance */
	//double *xa,*Pa;		/* fixed states and their covariance */
	netsol_t sol;		/* decoupled clock model solution */
	prcopt_t opt;		/* processing options */
	int nobs;			/* total obs number of all rcv */
	int fixedrec;		/* number of fixed receiver */
	int SSC[MAXSAT];	/* ambiguity datum for ssc */
	int SRC[MAXRCV];	/* ambiguity datum for src */
	int outc[MAXSAT];	/* outlock count */
	char time[32];		/* time string */
	int eflag;			/* epoch flag (0:first,1:mid,2:last) */

	gtime_t ctime;		/* extimated states time (GPST) */
	int nc,ncmax;		/* number of constraint states (with datum) */
	int nd,ndmax;		/* number of constraint states (without datum)*/
	csd_t *csd;			/* constraint data struct */
	dtm_t *dtm;			/* datum struct */
	double *DP;			/* constraint co-variance */

	//int skipfix;		/* epoch number of skipping fix (ar) */
} net_t;

typedef struct {        /* receiver parameter type */
	gtime_t time;
	rnxo_t rnx;			/* observation RINEX record */
	sta_t sta;			/* receiver information */

	int nobs;           /* obs number of this receiver in current epoch */
	obsd_t obs[MAXOBS];	/* observation data of this receiver in current epoch */

	int excrcv;         /* exclude this station without snx or blq */
	int datum[MAXSAT];	/* 0:not observed, 1:objected, 2:observed, 3:SSC, 4:SRC, 5:trans*/

	int ns;             /* number of valid satellites in thie station */
	double tt;          /* time difference between current and previous (s) */
	ssat_t ssat[MAXSAT];/* satellite status */

	double odisp[6*11]; /* ocean tide loading parameters */
	pcv_t pcvr;			/* receiver antenna parameters */
	int jumpc;

	int lock;			/* lock counter of rcv */
	unsigned int outc;	/* rcv outage counter */

	double pos[3];		/* blh coordinates after earth tides correction */
	double r[MAXOBS];	/* geometric distance and their variance */
	double rs[6*MAXOBS],varrs[MAXOBS];
	double dts[2*MAXOBS];
	double L[MAXOBS][NFREQ],Lc[MAXOBS];
	double P[MAXOBS][NFREQ],Pc[MAXOBS];
	double A4[MAXOBS];					/* MW combination */

	int fixed;			/* fixed rcv */
	ambt_t amb[MAXSAT];	/* ambiguity tracked struct */
	float ratio;		/* ratio value for AR */
	float rate;			/* success rate for AR */

	//int nofixcnt;		/* no fixed counter */
} rcv_t;

/* getMultiOBS.c -------------------------------------------------------------*/
extern pcvs_t *getpcvs_t();
extern obs_t *getobs_t();
extern nav_t *getnav_t();
extern rcv_t *getrcv_t();
extern int openMultiOBSFile(const char *file);
extern void getMultiOBSHeader(rcv_t *rcv, int n);
extern int getMultiOBSBody(gtime_t time, const char *opt, rcv_t *rcv, int n);
extern int getMultiOBS(rcv_t *rcv, const char *file);
extern void netRejSatAllRcv(rcv_t *rcv, int n, int rej_sat);
extern void netRejRcvAllSat(rcv_t *rcv);
extern void netRejRcvSat(const char *time, rcv_t *rcv, int n, int *outr, 
	int *outs);
extern void excludeObs(rcv_t *rcv,const prcopt_t *opt,const int n);
extern void outObsType(FILE *fp, const rcv_t *rcv, int n);
extern void outObsData(FILE *fp, const rcv_t *rcv, int n);

/* logDtm.c ------------------------------------------------------------------*/
extern int openDtmLog(const char* file);
extern void closeDtmLog(void);
extern void logDatumRcvPos(rcv_t *rcv,const int n);
extern void logDtm(gtime_t time,rcv_t *rcv,const int n);

/* netcmn.c ------------------------------------------------------------------*/
extern void initxNet(net_t *net, double xi, double var, int i);
extern int readSnxNet(const char *file, rcv_t *rcv, const int n);
extern void setPcvNet(gtime_t time, prcopt_t *popt, nav_t *nav, 
	const pcvs_t *pcvs, rcv_t *rcv, const int n);
extern void readOtlNet(const char *file, rcv_t *rcv, const int n);
extern int generateFinalClk(const prcopt_t *opt, const rcv_t *rcv, int n, 
	const char *tclkf, const char *atxf, const char *clkf);

/* setDatum.c ----------------------------------------------------------------*/
extern void setDatum(const net_t *uck, rcv_t *rcv, int n, int flag, int *SSC, 
	int *SRC);

#endif // !_NET_H
