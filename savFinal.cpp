//
// A YukawaPlasma class implementation
//   saves the final condition
//
// In-Gee Kim   August, 2015
//
// Separated from postEvolve(), Nov 25, 2015

// Private header
#include "YukawaPlasma.h"

void YukawaPlasma::savFinal(ofstream &logFile)
{
    //
    // YukawaPlasma class member function savFinal()
    //

    // Save final conditions
    cout << "Saving final conditions..." << endl;
    logFile << endl << "Saving final conditions..." << endl;

    //
    // save final velocities and poisions as iqp.out format
    // In-Gee Kim July, 2015
    //
    ofstream finalFile("initcb.out", ios::out);

    for (unsigned long i = 0; i < N_particles; i++) {
        if ( finalMode == 1 ) {
            // transform to the old format
            rx[i] += halfL;
            ry[i] += halfL;
            rz[i] += halfL;
        } // if ( finalMode == 1 )

        // save the information
        finalFile << scientific << setw(24) << setprecision(15)
            << vx[i] << " " << vy[i] << " " << vz[i] << " "
            << rx[i] << " " << ry[i] << " " << rz[i] << endl;
    } // for (unsigned long i = 0; i < N_particles; i++)

    finalFile.close();

    logFile << "Final velocities and positions are saved";
    if ( finalMode == 1 ) {
        logFile << " in old format" << endl;
    } else {
        logFile << " in new format" << endl;
    }

} // end of YukawaPlasma::savFinal()
