//
// An AutoCorrelations class implementation
//   save the averaged autocorrelation data into files
//
// Zv[typ][t] and Zj[t] contain the averaged velocity autocorrelation function
//  and the averaged current autocorrelation function, respectively.
//
// In-Gee Kim   December, 2015
//
// Separated from YukawaPlasma::postEvolve() Nov 25, 2015
// Separated from YukawaPlasma::calDiffusion() Dec 11, 2015
//

// Private headers
#include "AutoCorrelations.h"

void AutoCorrelations::savAutoCorr(void)
{
    // save the velocity auto-correlation function
    ofstream vacfFile("avg_vacf.out", ios::out); // for averaged vacf
    ofstream nvacFile("norm_vacf.out", ios::out); // for Zv(t)/Zv(0)
    for (unsigned long t = 0; t < M_steps; t++) {
        // bare vacf output
        vacfFile << right << fixed << setw(15) << setprecision(9)
            << ((double) t) * h << " "
            << scientific << setw(24) << setprecision(15);
        for (unsigned long typ = 0; typ < N_types; typ++) {
            vacfFile << Zv[typ][t] << " ";
        } // for (unsigned long typ = 0; typ < N_types; typ++)
        vacfFile << endl;

        // normalized vacf output
        nvacFile << right << fixed << setw(15) << setprecision(9)
            << ((double) t) * h << " "
            << scientific << setw(24) << setprecision(15);
        for (unsigned long typ = 0; typ < N_types; typ++) {
            nvacFile << (Zv[typ][t] / Zv[typ][0]) << " ";
        } // for (unsigned long typ = 0; typ < N_types; typ++)
        nvacFile << endl;
    } // for (unsigned long t = 0; t < M_steps; t++)

    vacfFile.close();
    nvacFile.close();

    // save the current auto-correlation function
    ofstream jacfFile("avg_jacf.out", ios::out);  // for averaged jacf
    ofstream njacFile("norm_jacf.out", ios::out); // for normalized jacf
    // ofstream nacfFile("norm_jacf.out", ios::out); // for Zj(t)/Zj(0)
    for (unsigned long t = 0; t < M_steps; t++) {
        // absolute jacf output
        jacfFile << right << fixed << setw(15) << setprecision(9)
        << ((double) t) * h << " "
        << scientific << setw(24) << setprecision(15);
        jacfFile << Zj[t] << endl;

        // normalized jacf output
        njacFile << right << fixed << setw(15) << setprecision(9)
        << ((double) t) * h << " "
        << scientific << setw(24) << setprecision(15);
        njacFile << (Zj[t] / Zj[0]) << endl;
    } // for (unsigned long t = 0; t < M_steps; t++)

    jacfFile.close();
    njacFile.close();

} // end of AutoCorrelations::savAutoCorr()
