
#include "changestats.users.h"

/* absdiffnodemix */

CHANGESTAT_FN(d_absdiffnodemix) {
  double change; Vertex t, h; int i, nnodes, nstats, statnum;
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    t = TAIL(i); h = HEAD(i);
    nnodes = INPUT_PARAM[0];
    nstats = INPUT_PARAM[1];
    change = fabs(INPUT_PARAM[t+1] - INPUT_PARAM[h+1]);
    for (statnum = 0; statnum < INPUT_PARAM[1]; statnum++)
      {
        if ((INPUT_PARAM[nnodes+t+1] == INPUT_PARAM[2*nnodes+statnum+2] &&
            INPUT_PARAM[nnodes+h+1] == INPUT_PARAM[2*nnodes+nstats+statnum+2]) ||
            (INPUT_PARAM[nnodes+t+1] == INPUT_PARAM[2*nnodes+nstats+statnum+2] &&
            INPUT_PARAM[nnodes+h+1] == INPUT_PARAM[2*nnodes+statnum+2]))
          {CHANGE_STAT[statnum] += IS_OUTEDGE(t,h) ? -change : change;}
      }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}


/* absdiffby */

CHANGESTAT_FN(d_absdiffby) {
  double change, offset, byval; Vertex t, h; int i;
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    t = TAIL(i); h = HEAD(i);
    byval = INPUT_PARAM[t + N_NODES];
    offset = INPUT_PARAM[0];
    if (byval == 1) {
      change = fabs(INPUT_PARAM[t] - INPUT_PARAM[h] - offset);
    } else {
      change = fabs(INPUT_PARAM[t] - INPUT_PARAM[h] + offset);
    }
    CHANGE_STAT[0] += IS_OUTEDGE(t,h) ? -change : change;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/* edist */

CHANGESTAT_FN(d_edist) {
  double change, d2, diff, expo, sign; Vertex t, h; int i, j, k, ndim, npow, coff;
  ndim = (int) INPUT_PARAM[0];
  npow = (int) INPUT_PARAM[1];
  coff = 2 + npow;  /* offset into INPUT_PARAM where coordinate blocks begin */
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    t = TAIL(i); h = HEAD(i);
    d2 = 0.0;
    for (k = 0; k < ndim; k++) {
      diff = INPUT_PARAM[coff + k * N_NODES + (t - 1)] - INPUT_PARAM[coff + k * N_NODES + (h - 1)];
      d2 += diff * diff;
    }
    sign = IS_OUTEDGE(t, h) ? -1.0 : 1.0;
    for (j = 0; j < npow; j++) {
      /* Euclidean distance raised to expo: (sqrt(d2))^expo == d2^(expo/2) */
      expo = INPUT_PARAM[2 + j];
      change = pow(d2, expo / 2.0);
      CHANGE_STAT[j] += sign * change;
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}


typedef struct {
  int *venues;
  int *lengths;
  int *positions;
  int binary;
} fuzzynodematch_storage;

I_CHANGESTAT_FN(i_fuzzynodematch) {
  ALLOC_STORAGE(1, fuzzynodematch_storage, sto);

  sto->venues = INTEGER(getListElement(mtp->R, "venues"));
  sto->lengths = INTEGER(getListElement(mtp->R, "lengths")) - 1;
  sto->positions = INTEGER(getListElement(mtp->R, "positions")) - 1;
  sto->binary = asInteger(getListElement(mtp->R, "binary"));
}

C_CHANGESTAT_FN(c_fuzzynodematch) {
  GET_STORAGE(fuzzynodematch_storage, sto);
  
  int i = 0;
  int j = 0;
  int *tailvenues = sto->venues + sto->positions[tail];
  int *headvenues = sto->venues + sto->positions[head];
  while(i < sto->lengths[tail] && j < sto->lengths[head]) {
    if(tailvenues[i] < headvenues[j]) {
      i++;
    } else if(tailvenues[i] > headvenues[j]) {
      j++;   
    } else {
      i++;
      j++;
      CHANGE_STAT[0]++;
    }
  }
  
  if(sto->binary) {
    CHANGE_STAT[0] = CHANGE_STAT[0] > 0;
  }
  
  if(edgestate) {
    CHANGE_STAT[0] = -CHANGE_STAT[0];
  }
}
