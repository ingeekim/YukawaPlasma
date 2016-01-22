//
// YukawaPlasma class implementation
//   post evolution properties
//
// In-Gee Kim   August, 2015
//

// Private header
#include "YukawaPlasma.h"

int YukawaPlasma::postEvolve(ofstream &logFile)
{
    //
    // YukawaPlasma class member function postEvolve()
    //
    // Invokes the molecular dynamics simulation indeed.
    //
    // It returns an int value 0 if the running has been successfully done,
    //   otherwise a certain error code, which will be defined later.
    //
    // Separated into a new function, Aug 14, 2015
    // Added the velocity auto-correlation function, Aug 17, 2015
    //
    // Transformed into a driver function for the post-processing routines.
    //  Nov 25, 2015


    // Save final conditions
    if ( finalMode ) {
        //
        // save final velocities and poisions as iqp.out format
        // In-Gee Kim July, 2015
        //
        savFinal(logFile);

        //
        // calculate and save the final velocity distributions
        //
        vel_dist(-2);

        logFile << "Final velocity distributions are saved" << endl;

        // save the radial distribution functions after normalization
        rad_dist(logFile);

        // calculate the diffusion coefficients by Green-Kubo theory
        calcDiffusion(logFile);
        if (idxEnsemble > 0) idxEnsemble++;

    } // if (finalMode)

    return idxEnsemble;
} // end of YukawaPlasma::postEvolve()
