//
// An AutoCorrelations class implementation
//   calculates the interdiffusion coefficients by the Fick's law
//
// In-Gee Kim   August, 2015
//
// Separated from YukawaPlasma::postEvolve() Nov 25, 2015
//

// Private header
#include "AutoCorrelations.h"

void AutoCorrelations::getFick(ofstream &logFile)
{
    //
    // Calculate the inter-diffusion coefficient D^0_12
    //  in units of a_1^2 \omega_1, where D^0_12 = \int Z_j(t) / 3 x_l x_h
    //
    // In-Gee Kim, December 2015

    double inter_Diff = 0.0;
    // End points
    inter_Diff = (half * (Zj[0] + Zj[M_steps - 1]));

    // Mid-points
    for (unsigned long t = 1; t < (M_steps - 1); t++) {
        inter_Diff += Zj[t];
    } // for (unsigned long t = 1; t < (M_steps - 1); t++)

    // integral
    inter_Diff *= h;

    // multiply the final factdors for the unit
    inter_Diff *= oneThird;
    inter_Diff /= ((double) N_ptls);
    inter_Diff /= (x_l * x_h);

    // inter-diffusion coefficient from the current auto-correlation
    logFile << endl;
    logFile << "The interdiffusion coefficient "
        <<  "from the current auto-correlation function" << endl;
    logFile << "\tD^0_ij = "
        << inter_Diff << " / a_l^2 w_l = "
        << inter_Diff * a_l * a_l * w_l << " cm^2 / s" << endl;

} // end of AutoCorrelations::getFick()
