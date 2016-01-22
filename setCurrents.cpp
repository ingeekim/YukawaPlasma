//
// An AutoCorrelations class implementation
//   calculates the current density
//
// In-Gee Kim   December, 2015
//
// Separated from postEvolve() Nov 25, 2015
// Separated from calcDiffusion() Dec 11, 2015

// Private header
#include "AutoCorrelations.h"

void AutoCorrelations::setCurrents(void)
{
    // calculating current density
    for (unsigned long t = 0; t < M_steps; t++) {
        for (unsigned long i = 0; i < N_light; i++) {
            jxt[t] += x_l * vxt[t][i];
            jyt[t] += x_l * vyt[t][i];
            jzt[t] += x_l * vzt[t][i];
        } // for (unsigned long i = 0; i < N_light; i++)

        for (unsigned long i = 0; i < N_heavy; i++) {
            jxt[t] -= x_h * vxt[t][i + N_light];
            jyt[t] -= x_h * vyt[t][i + N_light];
            jzt[t] -= x_h * vzt[t][i + N_light];
        } // for (unsigned long i = 0; i < N_heavy; i++)
    } // for (unsigned long t = 0; t < (M_timeStep / M_snapShot); t++)

} // end of AutoCorrelations::calcCurrents()
