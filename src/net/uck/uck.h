
#ifndef _UCK_H
#define _UCK_H

#include "net/net.h"

#define UNF(opt)   ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)	 /* number of frequency */
#define UNSC(opt)  (MAXSAT)										 /* number of satellite clock */
#define UNSPB(opt) ((opt)->ionoopt==IONOOPT_IFLC?0:UNF(opt)*MAXSAT)/* number of satellite phase bias */
#define UNSCB(opt) ((opt)->scb?UNF(opt)*MAXSAT:0)				   /* number of satellite code bias */
#define UNSAT(opt) (UNSC(opt)+UNSPB(opt)+UNSCB(opt))

#define UNRC(opt)   (1)											/* number of receiver clock for single receiver */
#define UNRCB(opt)  ((opt)->rcb?UNF(opt):0)						/* number of receiver code bias drift */
#define UNT(opt)	((opt)->tropopt==TROPOPT_EST?1:3)			/* number of troposphere delay for single receiver */
#define UNI(opt)	((opt)->ionoopt==IONOOPT_EST?MAXSAT:0)		/* number of ionosphere delay for single receiver */
#define UNAMB(opt)	(UNF(opt)*MAXSAT)							/* number of ambiguity for single receiver */
#define UNRPB(opt)  ((opt)->ionoopt==IONOOPT_IFLC?0:UNF(opt))	/* number of receiver phase bias for single receiver */
#define UNR(opt)	(UNRC(opt)+UNRCB(opt)+UNT(opt)+UNI(opt)+UNAMB(opt)+UNRPB(opt))	/* number of parameters for single receiver */

/* s for satellite number (from 1), r for receiver number (from 1), f for frequency number (from 1) */
#define UISC(s,opt)		(s-1)										/* index of satellite clock */
#define UISPB(s,f,opt)	(UNSC(opt)+(s-1)*UNF(opt)+(f-1))			/* index of satellite phase bias */
#define UISCB(s,f,opt)	(UNSC(opt)+UNSPB(opt)+(s-1)*UNF(opt)+(f-1))	/* index of satellite code bias */

#define UIRC(r,opt)		 (UNSAT(opt)+(r-1)*UNR(opt))						/* index of receiver clock */
#define UIRCB(r,f,opt)	 (UNSAT(opt)+(r-1)*UNR(opt)+UNRC(opt)+(f-1))		/* index of receiver code bias */
#define UIT(r,opt)		 (UNSAT(opt)+(r-1)*UNR(opt)+UNRC(opt)+UNRCB(opt))	/* index of troposphere delay */
#define UII(r,s,opt)	 (UNSAT(opt)+(r-1)*UNR(opt)+UNRC(opt)+UNRCB(opt)+UNT(opt)+(s-1))			/* index of ionosphere delay */
#define UIAMB(r,s,f,opt) (UNSAT(opt)+(r-1)*UNR(opt)+UNRC(opt)+UNRCB(opt)+UNT(opt)+UNI(opt)+UNF(opt)*(s-1)+(f-1))/* index of ambiguity */
#define UIRPB(r,f,opt)	 (UNSAT(opt)+(r-1)*UNR(opt)+UNRC(opt)+UNRCB(opt)+UNT(opt)+UNI(opt)+UNAMB(opt)+(f-1))	/* index of receiver phase bias */



//typedef struct {
//	gtime_t time;		/* time (GPST) */
//	int stat;
//	float ratio;        /* AR ratio factor for valiation */
//} ucksol_t;

//typedef struct {
//	int nx;				/* number of float states */
//	double *x,*P;		/* float states and their covariance */
//	ucksol_t sol;		/* decoupled clock model solution */
//	prcopt_t opt;		/* processing options */
//	int nobs;			/* total obs number of all rcv */
//	int fixedrec;		/* number of fixed receiver */
//	int SSC[MAXSAT];	/* ambiguity datum for ssc */
//	int SRC[MAXRCV];	/* ambiguity datum for src */
//	int outc[MAXSAT];	/* sat out counter */
//	char time[32];		/* time string */
//	int eflag;			/* epoch flag (0:first,1:mid,2:last) */
//	
//	gtime_t ctime;		/* extimated states time (GPST) */
//	int nc,ncmax;		/* number of constraint states (with datum) */
//	int nd,ndmax;		/* number of constraint states (without datum)*/
//	csd_t *csd;			/* constraint data struct */
//	dtm_t *dtm;			/* datum struct */
//	double *DP;			/* constraint co-variance */
//} net_t;

//typedef struct {        /* receiver parameter type */
//	gtime_t time;
//	rnxo_t rnx;			/* observation RINEX record */
//	sta_t sta;			/* receiver information */
//
//	int nobs;           /* obs number of this receiver in current epoch */
//	obsd_t obs[MAXOBS];	/* observation data of this receiver in current epoch */
//
//	int excrcv;         /* exclude this station without snx or blq */
//	int datum[MAXSAT];	/* 0:not observed, 1:objected, 2:observed, 3:SSC, 4:SRC, 5:trans*/
//	
//    int ns;             /* number of valid satellites in thie station */
//    double tt;          /* time difference between current and previous (s) */
//    ssat_t ssat[MAXSAT];/* satellite status */
//
//	double odisp[6*11]; /* ocean tide loading parameters */
//	pcv_t pcvr;			/* receiver antenna parameters */
//	int jumpc;
//
//	int lock;			/* lock counter of rcv */
//	unsigned int outc;	/* rcv outage counter */
//
//	double pos[3];		/* blh coordinates after earth tides correction */
//	double r[MAXOBS];	/* geometric distance and their variance */
//    double rs[6*MAXOBS],varrs[MAXOBS];
//	double dts[2*MAXOBS];
//	double L[MAXOBS][NFREQ],Lc[MAXOBS];
//	double P[MAXOBS][NFREQ],Pc[MAXOBS];
//
//	int fixed;			/* fixed rcv */
//	ambt_t amb[MAXSAT];	/* ambiguity tracked struct */
//	float ratio;		/* ratio value for AR */
//	float rate;			/* success rate for AR */
//} rcv_t;

/* functions -----------------------------------------------------------------*/

/* uckPost.c -----------------------------------------------------------------*/
extern int uckPost(const prcopt_t *popt, const solopt_t *sopt, 
	const filopt_t *fopt);
/* outProduct_u.c ------------------------------------------------------------*/
extern void outProduct_u(FILE *fp, const net_t *net, const rcv_t *rcv, int n);

/* estProduct_u.c ------------------------------------------------------------*/
extern int estProduct_u(net_t *uck, rcv_t *rcv, int n, const nav_t *nav, 
	int flag);
/* udStates_u.c --------------------------------------------------------------*/
extern void udStates_u(net_t *uck, rcv_t *rcv, int n, const nav_t *nav);
/* setDatum_u.c --------------------------------------------------------------*/
//extern void setDatum_u(const net_t *uck, rcv_t *rcv, int n, int flag, int *SSC, 
//	int *SRC);
/* uckRes.c ------------------------------------------------------------------*/
extern int uckRes(int post, const net_t *uck, rcv_t *rcv, int n, 
	const nav_t *nav, const double *x, const double *P, int *ix, double *vl, 
	double *HL, double *DL, double *vp, double *HP, double *DP, int *nc, 
	int flag);
/* uckAmbResol.c -------------------------------------------------------------*/
extern void udTrackAmb_u(net_t *uck, rcv_t *rcv, int n, int ar);
extern void udAmbDatum_u(const prcopt_t *opt, rcv_t *rcv, int n);
extern int uckAmbResol(const net_t *uck, rcv_t *rcv, int n, double *x, 
	double *P);

/* logRes_u.c ----------------------------------------------------------------*/
extern int openResLog_u(const char *file);
extern void closeResLog_u(void);
extern void logRes_u(net_t *uck,rcv_t *rcv,const int n);
/* logZTD.c ------------------------------------------------------------------*/
extern int openZtdLog_u(const char *file);
extern void closeZtdLog_u(void);
extern void logZtd_u(net_t *uck,rcv_t *rcv,const int n);
/* logION.c ------------------------------------------------------------------*/
extern int openIonLog_u(const char *file);
extern void closeIonLog_u(void);
extern void logIon_u(net_t *uck,rcv_t *rcv,const int n);
/* logAmb_u.c ----------------------------------------------------------------*/
extern int openAmbLog_u(const char *file);
extern void closeAmbLog_u(void);
extern void logAmb_u(net_t *uck,rcv_t *rcv,const int n);
/* logPhb.c ------------------------------------------------------------------*/
extern int openUPDLog_u(const char *file);
extern void closeUPDLog_c(void);
extern void logUPD_u(net_t *uck, rcv_t *rcv, const int n);
/* logCob_u.c ----------------------------------------------------------------*/
extern int openCobLog_u(const char *file);
extern void closeCobLog_u(void);
extern void logCob_u(net_t *uck, rcv_t *rcv, const int n);

/* record_u.c ----------------------------------------------------------------*/
extern void recordOutFile_u(const char *file);
extern void openRecord_u(int lastEp);
extern void closeRecord_u(int lastEp);
extern void record_u(int type, const char *time, const char *rcv, int sat, 
	int frq, int trp, int datum, double value, int lastEp);
extern void recordVar_u(const double *DP, int nc, int lastEp);
extern int getRecord_u(net_t *uck, const char *fileIndex, const char *fileValue);
extern int copyConstraint_u(net_t *uck, rcv_t *rcv, int n);

#endif /* _UCK_H */