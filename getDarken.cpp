//
// An AutoCorrelations class implementation
//   calculates the diffusion coefficients by Darken rule
// The results are printed into the logFile
//
// In-Gee Kim   December, 2015
//
// Separated from YukawaPlasma::postEvolve() Nov 25, 2015
// Separated from YukawaPlasam::calcDiffusion() Dec 11, 2015
//

// Private header
#include "AutoCorrelations.h"

void AutoCorrelations::getDarken(ofstream &logFile)
{
    //
    // Calculate the self-diffusion coefficienets D[i] (i = 0, ..., N_type)
    //   in units of a_1^2 \omega_1, where D = \int Z(t) dt / 3
    // The time integration is performed over the finite data set Z[i][t]
    //   by using the trapezoidal integration
    //
    // In-Gee Kim, August 2015

    // define the self-diffusion coefficients
    vector<double> self_Diff(N_types, 0.0);

    for (unsigned long typ = 0; typ < N_types; typ++) {
        // End points
        self_Diff[typ] = (half * (Zv[typ][0] + Zv[typ][M_steps - 1]));

        // Mid-points
        for (unsigned long t = 1; t < (M_steps - 1); t++) {
            self_Diff[typ] += Zv[typ][t];
        }

        // integral
        self_Diff[typ] *= h;

        // Final factor from the unit
        self_Diff[typ] *= oneThird;
    } // for (unsigned long typ = 0; typ < N_types; typ++)

    // record the self-diffusion coefficients to logFile
    logFile << endl << "The diffusion coefficients are" << endl;
    for (unsigned long typ = 0; typ < N_types; typ++) {
        logFile << "\tDiff[" << typ << "] = " << self_Diff[typ]
        << " / a_l^2  w_l = "
        << self_Diff[typ] * a_l * a_l * w_l << " cm^2 / s" << endl;
    } // for (unsigned long typ = 0; typ < N_types; typ++)
    logFile << endl;

    //
    // inter-diffusion coefficient by Darken rule
    //
    //  In-Gee Kim, November 2015
    double diff_Darken = 0.0;
    diff_Darken = (x_h * self_Diff[0] + x_l * self_Diff[1]);
    logFile << "The interdiffusion coefficient by Darken rule" << endl;
    logFile << "\tD^0_ij = "
        << diff_Darken << " / a_l^2 w_l = "
        << diff_Darken * a_l * a_l * w_l << " cm^2 / s" << endl;

} // end of AutoCorrelations::getDarken()
