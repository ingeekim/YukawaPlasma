//
// A YukawaPlasma class implementation
//   calculates the diffusion coefficients
//
// In-Gee Kim   August, 2015
//
// Separated from postEvolve() Nov 25, 2015
//

// Private headers
#include "YukawaPlasma.h"
#include "AutoCorrelations.h"

void YukawaPlasma::calcDiffusion(ofstream &logFile)
{
    //
    // calculate and save the velocity auto-correlation functions
    //


    // declare an object of AutoCorrelations
    AutoCorrelations autoCorr(N_samples, idxEnsemble,
                              M_timeSteps, M_snapShots, N_particles);
    autoCorr.setSpecies(N_light, N_heavy, x_h, x_l);
    double t_s = (dt * ((double) M_snapShots));
    autoCorr.setUnits(t_s, a_l, w_l);

    // read the velocity data
    autoCorr.readVelocities();

    // calculating the current density
    autoCorr.setCurrents();
    logFile << "Calculating velocity auto-correlation functions..." << endl;

    // calculating the autocorrelation functions
    autoCorr.calcAutoCorr();
    logFile << "Autocorrelation functions are calculated." << endl;

    // accumulate the autocorrelation functions
    autoCorr.accAutoCorr();
    logFile << "Autocorrelation functions are accumulated." << endl;

    // average the autocorrelation functions
    autoCorr.avgAutoCorr();
    logFile << "Autocorrelation functions are averaged." << endl;

    // save the averaged autocorrelation functions
    autoCorr.savAutoCorr();
    logFile << "The velocity auto-correlation functions are saved." << endl;

    // check the statistical fitness
    int isFit = 0;
    isFit = autoCorr.chkStatFit(logFile);

    // calculate the diffusion coefficients in Darken rule
    autoCorr.getDarken(logFile);

    // calculate the inter-diffusion coefficient D^0_12
    autoCorr.getFick(logFile);

} // end of YukawaPlasma::calcDiffusion()
