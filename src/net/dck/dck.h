
#ifndef _DCK_H
#define _DCK_H

#include "net/net.h"

#define DNF(opt)	 2
#define CLK_NUM(opt) 3
#define DNSC(opt)	(CLK_NUM(opt)*MAXSAT)				/* number of satellite clock */
#define DNRC(opt)	(CLK_NUM(opt))						/* number of receiver clock for single station */
#define DNT(opt)	((opt)->tropopt==TROPOPT_EST?1:3)	/* number of troposphere delay for single station */
#define DNAMB(opt)	(DNF(opt)*MAXSAT)					/* number of ambiguity for single station */
#define DNR(opt)	(DNRC(opt)+DNT(opt)+DNAMB(opt))		/* number of parameters except satellite clocks for single station */

/* s for satellite number (from 1), r for station number (from 1), n for equation number (from 1) */
#define DISC(s,n,opt)		((s-1)*CLK_NUM(opt)+(n-1))				/* index of satellite clock or clock like */
#define DIRC(r,n,opt)		(DNSC(opt)+(r-1)*DNR(opt)+(n-1))			/* index of receiver clock or clock like */
#define DIT(r,opt)			(DNSC(opt)+(r-1)*DNR(opt)+DNRC(opt))	/* index of tropsphere delay */
#define DIAMB(r,s,n,opt)	(DNSC(opt)+(r-1)*DNR(opt)+DNRC(opt)+DNT(opt)+DNF(opt)*(s-1)+(n-1))	/* index of ambiguity */

/* functions -----------------------------------------------------------------*/

/* dckPost.c -----------------------------------------------------------------*/
extern int dckPost(const prcopt_t *popt, const solopt_t *sopt, 
	const filopt_t *fopt);

/* outProduct_d.c ------------------------------------------------------------*/
extern void outProduct_d(FILE *fp, const net_t *dck, const rcv_t *rcv, int n);

/* estProduct_d.c ------------------------------------------------------------*/
extern int estProduct_d(net_t *dck, rcv_t *rcv, int n, const nav_t *nav, 
	int flag);

/* udStates_d.c --------------------------------------------------------------*/
extern void udStates_d(net_t *dck, rcv_t *rcv, int n);

/* dckRes.c ------------------------------------------------------------------*/
extern int dckRes(int post, const net_t *dck, rcv_t *rcv, int n,
	const nav_t *nav, const double *x, const double *P, int *ix, double *vl, 
	double *HL, double *DL, double *vp, double *HP, double *DP, int *nc, 
	int tnv, int flag, int *outr, int *outs);

/* dckNetAr.c ----------------------------------------------------------------*/
extern void udAmbDatum_d(const prcopt_t *opt, rcv_t *rcv, int n);
extern void udTrackAmb_d(net_t *dck, rcv_t *rcv, int n, int ar);
extern int dckNetAr(const net_t *dck, rcv_t *rcv, int n, double *x, double *P);

/* logRes_d.c ----------------------------------------------------------------*/
extern int openResLog_d(const char *file);
extern void closeResLog_d(void);
extern int logRes_d(net_t *dck,rcv_t *rcv,const int n);

/* logZtd_d.c ----------------------------------------------------------------*/
extern int openZtdLog_d(const char *file, int opt);
extern void closeZtdLog_d(void);
extern void logZtd_d(net_t *dck,rcv_t *rcv,const int n);

/* logAmb_d.c ----------------------------------------------------------------*/
extern int openAmbLog_d(const char *file);
extern void closeAmbLog_d(void);
extern void logAmb_d(net_t *dck,rcv_t *rcv,const int n);

/* record_d.c ----------------------------------------------------------------*/
extern void recordOutFile_d(const char *file);
extern void openRecord_d(int lastEp);
extern void closeRecord_d(int lastEp);
extern void record_d(int type, const char *time, const char *rcv, int sat, 
	int frq, int trp, int datum, double value, int lastEp);
extern void recordVar_d(const double *DP, int nc, int lastEp);
extern int getRecord_d(net_t *uck, const char *fileIndex, const char *fileValue);
extern int copyConstraint_d(net_t *uck, rcv_t *rcv, int n);

#endif /* _DCK_H */