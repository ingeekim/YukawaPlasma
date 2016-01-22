//
// An AutoCorrelations class implementation
//   check the statistical fitness of the averaged autocorrelation functions
// When the criteria is satisfied,
//    the function set the idxEnsemble to be negative value and returns it.
//
// In-Gee Kim   December, 2015
//

// System header
#include <cmath>
using std::sqrt;

// Private header
#include "AutoCorrelations.h"

int AutoCorrelations::chkStatFit(ofstream &logFile)
{
    // avoid the statistical calculations with single sample
    if (idxEnsemble == 1) return idxEnsemble;

    // open the history files
    ifstream hist_vacfFile("hist_vacf.out", ios::in);
    ifstream hist_jacfFile("hist_jacf.out", ios::in);

    // load data over the samples in the ensemble
    double time = 0.0;
    vector<double> sigmaV0(M_steps, 0.0);
    vector<double> sigmaV1(M_steps, 0.0);
    vector<double> sigmaJ(M_steps, 0.0);
    double iZv0 = 0.0;
    double iZv1 = 0.0;
    double iZj = 0.0;
    for (int id = 1; id <= idxEnsemble; id++) {
        // loop over time of the sample
        for(unsigned long t = 0; t < M_steps; t++) {
            // load data
            hist_vacfFile >> time >> iZv0 >> iZv1;
            hist_jacfFile >> time >> iZj;

            // sum the square distances from the averaged ones over the samples
            sigmaV0[t] += (iZv0 - Zv[0][t]) * (iZv0 - Zv[0][t]);
            sigmaV1[t] += (iZv1 - Zv[1][t]) * (iZv1 - Zv[1][t]);
            sigmaJ[t] += (iZj - Zj[t]) * (iZj - Zj[t]);
        } // for(unsigned long t = 0; t < M_steps; t++)
    } // for (int id = 1; id <= idxEnsemble; id++)

    // release the history files
    hist_vacfFile.close();
    hist_jacfFile.close();

    // obtain the sigma(t)
    double denom = 1.0;
    denom = ((double) idxEnsemble);
    for (unsigned long t = 0; t < M_steps; t++) {
        // average the mean-squares
        sigmaV0[t] /= denom;
        sigmaV1[t] /= denom;
        sigmaJ[t] /= denom;

        // root-mean-square
        sigmaV0[t] = sqrt(sigmaV0[t]);
        sigmaV1[t] = sqrt(sigmaV0[t]);
        sigmaJ[t] = sqrt(sigmaJ[t]);
    } // for (unsigned long t = 0; t < M_steps; t++)

    // find the maximum standard deviations
    double max_sigmaV0 = 0.0;
    double max_sigmaV1 = 0.0;
    double max_sigmaJ = 0.0;
    for (unsigned long t = 0; t < M_steps; t++) {
        if (sigmaV0[t] > max_sigmaV0) max_sigmaV0 = sigmaV0[t];
        if (sigmaV1[t] > max_sigmaV1) max_sigmaV1 = sigmaV1[t];
        if (sigmaJ[t] > max_sigmaJ) max_sigmaJ = sigmaJ[t];
    } // for (unsigned long t = 0; t < M_steps; t++)

    // print the maximum standard deviations
    logFile << endl;
    logFile << "The maximum standard deviation for \n";
    logFile << "\tvacf[0] = " << max_sigmaV0 << " (a_l w_l)^2" << endl;
    logFile << "\tvacf[1] = " << max_sigmaV1 << " (a_l w_l)^2" << endl;
    logFile << "\t   jacf = " << max_sigmaJ << " (a_l w_l)^2." << endl;
    logFile << endl;

    // check the statistical fits
    const double criterion = 0.4;
    int flag = 0;
    if (max_sigmaV0 < criterion) flag--;
    if (max_sigmaV1 < criterion) flag--;
    if (max_sigmaJ < criterion) flag--;

    if (flag == -3) idxEnsemble = flag;

    return idxEnsemble;
} // end of AutoCorrelations::calcCurrents()
