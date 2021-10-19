#ifndef _QC_H
#define _QC_H

#include "base/base.h"

extern void detClockJump(ssat_t *ssat, int *jumpc, obsd_t* obs, const int n, 
    const nav_t* nav, const char *rcv);
extern void detCycleSlip(const prcopt_t *opt, ssat_t *ssat, obsd_t *obs, 
    const int n, const nav_t *nav, const char *rcv, const double tt, int *jumpc);


#endif /* _QC_H */
