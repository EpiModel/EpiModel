
#ifndef CHANGESTATS_H
#define CHANGESTATS_H

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_storage.h"
#include "ergm_Rutil.h"

CHANGESTAT_FN(d_absdiffnodemix);
CHANGESTAT_FN(d_absdiffby);
CHANGESTAT_FN(d_edist);
C_CHANGESTAT_FN(c_fuzzynodematch);

#endif
