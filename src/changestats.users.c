
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
