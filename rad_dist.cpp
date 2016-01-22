//
// A YukawaPlasma class implementation
//   calculates the radial distributio functions
//
// In-Gee Kim   August, 2015
//
// Separated from postEvolve(), Nov 25, 2015

// Private header
#include "YukawaPlasma.h"

void YukawaPlasma::rad_dist(ofstream &logFile)
{
    //
    // YukawaPlasma class member function rad_dist()
    //

    const unsigned long M_steps = (M_timeSteps / M_snapShots);
    
    // save the radial distribution functions after normalization
    ofstream rdfFile("rdf.out", ios::out);

    double x_ll = ((double) N_light) / (edgeL * edgeL * edgeL);
    x_ll *= dr * ((double) N_light) * ((double) M_steps);
    double x_hh = ((double) N_heavy) / (edgeL * edgeL * edgeL);
    x_hh *= dr * ((double) N_heavy) * ((double) M_steps);
    double x_lh = ((double) N_light) / (edgeL * edgeL * edgeL);
    x_lh *= dr * ((double) N_heavy) * ((double) M_steps);
    double r = 0.0;
    double v0 = 0.0;
    for (unsigned long i = 0; i < R_bins; i++) {
        r = (((double) i) + half) * dr;
        v0 = fourPi * (r * r + r * dr + oneThird * dr * dr);
        g_ll[i] /= (v0 * x_ll);
        g_hh[i] /= (v0 * x_hh);
        g_lh[i] /= (2.0 * v0 * x_lh);

        rdfFile << right << fixed << setw(15) << setprecision(9)
                << r << " "
                << scientific << setw(24) << setprecision(15)
                << g_ll[i] << " " << g_hh[i] << " " << g_lh[i] << endl;
    }

    rdfFile.close();

    logFile << "The radial distribution functions are saved" << endl;

} // end of YukawaPlasma::postEvolve()
