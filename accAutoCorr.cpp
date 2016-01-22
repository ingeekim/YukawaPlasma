//
// An AutoCorrelations class implementation
//   accumulates the autocorrelation functions
//   and save the functions into the corresponding files
//
// In-Gee Kim   December, 2015
//

// Private header
#include "AutoCorrelations.h"

void AutoCorrelations::accAutoCorr(void)
{
    //
    // prepare the files and data
    //
    if (idxEnsemble > 0) {
        cerr << "Load the existing accumulated autocorrelation functions." << endl;
        // open the accumulated autocorrelation files
        ifstream acvacfFile("acc_vacf.out", ios::in);
        ifstream acjacfFile("acc_jacf.out", ios::in);

        // load the accumulated autocorrelation functions from the files
        double time = 0.0;
        for (unsigned long t = 0; t < M_steps; t++) {
            acvacfFile >> time;
            acjacfFile >> time;
            for (unsigned long ty = 0; ty < N_types; ty++) {
                acvacfFile >> accZv[ty][t];
            } // for (unsigned long ty = 0; ty < N_types; ty++)
            acjacfFile >> accZj[t];
        } // for (unsigned long t = 0; ty < M_steps; t++)

        // release the ifstremed accumulated autocorrelation files
        acvacfFile.close();
        acjacfFile.close();
    }

    //
    // now create the accumulated autocorrelation files for saving
    //
    ofstream acc_vacfFile("acc_vacf.out", ios::out);
    ofstream acc_jacfFile("acc_jacf.out", ios::out);

    // accumulate the accumulated autocorrelation functions
    for (unsigned long t = 0; t < M_steps; t++) {
        for (unsigned long ty = 0; ty < N_types; ty++) {
            accZv[ty][t] += Zv[ty][t];
        } // for (unsigned long ty = 0; t < N_types; ty++)
        accZj[t] += Zj[t];
    } // for (unsigned long t = 0; t < M_steps; t++)

    // save the accumulated autocorrelation functions
    double time = 0.0;
    for (unsigned long t = 0; t < M_steps; t++) {
        time = h * ((double) t);
        acc_vacfFile << right << fixed << setw(15) << setprecision(9)
            << time << " "
            << scientific << setw(24) << setprecision(15);
        for (unsigned long ty = 0; ty < N_types; ty++) {
            acc_vacfFile << accZv[ty][t] << " ";
        } // for (unsigned long ty = 0; ty < N_types; ty++)
        acc_vacfFile << endl;

        acc_jacfFile << right << fixed << setw(15) << setprecision(9)
            << time << " "
            << scientific << setw(24) << setprecision(15);
        acc_jacfFile << accZj[t] << endl;
    } // for (unsigned long t = 0; t < M_steps; t++)

    // release the accumulated autocorrelation function files
    acc_vacfFile.close();
    acc_jacfFile.close();

    //
    // create the history files
    //
    ofstream hist_vacfFile("hist_vacf.out", ios::app);
    ofstream hist_jacfFile("hist_jacf.out", ios::app);

    // save the autocorrelation functions to the histofy files
    time = 0.0;
    for (unsigned long t = 0; t < M_steps; t++) {
        time = h * ((double) t);
        hist_vacfFile << right << fixed << setw(15) << setprecision(9)
            << time << " "
            << scientific << setw(24) << setprecision(15);
        for (unsigned long ty = 0; ty < N_types; ty++) {
            hist_vacfFile << Zv[ty][t] << " ";
        } // for (unsigned long ty = 0; ty < N_types; ty++)
        hist_vacfFile << endl;

        hist_jacfFile << right << fixed << setw(15) << setprecision(9)
            << time << " "
            << scientific << setw(24) << setprecision(15);
        hist_jacfFile << accZj[t] << endl;
    } // for (unsigned long t = 0; t < M_steps; t++)

    // release the history files
    hist_vacfFile.close();
    hist_jacfFile.close();


} // end of AutoCorrelations::calcCurrents()
