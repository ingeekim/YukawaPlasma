//
// An AutoCorrelations class implementation
//   average the autocorrelation functions
//
// In-Gee Kim   December, 2015
//

// Private header
#include "AutoCorrelations.h"

void AutoCorrelations::avgAutoCorr(void)
{

    // average the autocorrelation functions
    double denom = 1.0;
    denom = ((double) idxEnsemble);
    for (unsigned long t = 0; t < M_steps; t++) {
        for (unsigned long ty = 0; ty < N_types; ty++) {
            Zv[ty][t] = (accZv[ty][t] / denom);
        } // for (unsigned long ty = 0; t < N_types; ty++)
        Zj[t] = (accZj[t] / denom);
    } // for (unsigned long t = 0; t < M_steps; t++)

} // end of AutoCorrelations::calcCurrents()
