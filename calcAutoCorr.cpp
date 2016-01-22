//
// An AutoCorrelations class implementation
//   calculates the autocorrelation functions
//
// In-Gee Kim   December, 2015
//
// Separated from postEvolve() Nov 25, 2015
// Separated from calcDiffusion() Dec 11, 2015
//

// Private header
#include "AutoCorrelations.h"

void AutoCorrelations::calcAutoCorr(void)
{

    long counter = -1;
    double partialSum = 0.0;
    double scalarProduct = 0.0;
    for (unsigned long i = 0; i < M_steps; i++) {
        // loop for the different time origin
        counter = -1; // reset the counter for each origin
        for (unsigned long j = i; j < M_steps; j++) {
            // loop for the given time series
            counter++;

            // for the velocity auto-correlation functions
            partialSum = 0.0;
            for (unsigned long k = 0; k < N_light; k++) {
                // light species
                scalarProduct = vxt[j][k] * vxt[i][k]
                    + vyt[j][k] * vyt[i][k] + vzt[j][k] * vzt[i][k];
                partialSum += scalarProduct;
            } // for (unsigned long k = 0; k < N_light; k++)
            Zv[0][counter] += (partialSum / ((double) N_light));

            partialSum = 0.0;
            for (unsigned long k = N_light; k < N_ptls; k++) {
                // heavy species
                scalarProduct = vxt[j][k] * vxt[i][k]
                    + vyt[j][k] * vyt[i][k] + vzt[j][k] * vzt[i][k];
                partialSum += scalarProduct;
            } // for (unsigned long k = N_light; k < N_ptls; k++)
            Zv[1][counter] += (partialSum / ((double) N_heavy));

            // for the current auto-correlation function
            partialSum = 0.0;
            scalarProduct = jxt[j] * jxt[i] + jyt[j] * jyt[i] + jzt[j] * jzt[i];
            Zj[counter] += scalarProduct;
        } // for (unsigned long j = 0; j < M_steps; j++)
    } // for (unsigned long i = 0; i < M_steps; i++)

    // normalize velocity auto-correlation functions
    double denominator = ((double) (M_steps - 1));
    for (unsigned long typ = 0; typ < N_types; typ++) {
        for (unsigned long t = 0; t < M_steps; t++) {
            Zv[typ][t] /= denominator; // time average
        } // for (unsigned long i = 0; i < M_steps; i++)
    } // for (unsigned long typ = 0; typ < N_types; typ++)

    // normalize the current auto-correlation functions
    for (unsigned long t = 0; t < M_steps; t++) {
        Zj[t] /= denominator;
    } // for (unsigned long t = 0; t < M_steps; t++


} // end of AutoCorrelations::postEvolve()
