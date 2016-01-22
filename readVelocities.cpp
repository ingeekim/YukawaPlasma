//
// An AutoCorrelations class implementation
//   read the velocities from the file, "iqp.out".
//
// In-Gee Kim   December, 2015
//
// Separated from postEvolve() Nov 25, 2015
// Separated from YukawaPlasa::calDiffusion() Dec 11, 2015
//

// Private header
#include "AutoCorrelations.h"

void AutoCorrelations::readVelocities(void)
{
    // read data from the file 'iqp.out'
    ifstream qpFile("iqp.out", ios::in);

    double discard = 0.0;
    for (unsigned long t = 0; t < M_steps; t++) {
        for (unsigned long i = 0; i < N_ptls; i++) {
            qpFile >> vxt[t][i] >> vyt[t][i] >> vzt[t][i]
                    >> discard >> discard >> discard
                    >> discard >> discard >> discard;
        } // for (unsigned long i = 0; i < N_particles; i++)
    } // for (unsigned long t = 0; t < (M_timeStep / M_snapShot); t++)

    qpFile.close();

}
