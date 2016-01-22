//
// The Yukawa binary ionic mixture molecular dynamics simulation program
//
// In-Gee Kim, July, 2015
//

//
// Declare the system headers
//
#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

#include <fstream>
using std::ifstream;

#include <cstdlib>
using std::exit;

#include <string>
using std::string;

//
// Declare the user headers
# include "YukawaPlasma.h"
//

//
// Name spaces
//
using namespace std;

int main(int argc, char *argv[])
{
    cout << "#####\n"
        << "ybim: Yukawa binary ionic mixtures molecular dynamics\n"
        << endl;

    int N_samples = 1;
    if (argc > 1) {  // check the arguments
        // Initiate the time decomposition method
        // usage: ybim 123
        //    The number 123 is the maximum number of samples
        //    for ensemble average the autocorrelation functions.
        // In-Gee Kim, December, 2015

        cout << "The time decomposition method is invoked." << endl;
        string str_samples = argv[--argc];
        N_samples = stoi(str_samples);
        cout << "\tThe maximum number of samples in the ensemble are "
            << N_samples << endl;
    } // end of the arguments check

    // Declare the ybim object
    // The input processing and initialization will be done.
    ofstream logFile("diagn.out", ios::out);
    YukawaPlasma ybim(logFile);
    ybim.set_ensemble(N_samples);

    // Let ybim evolve in time.
    int status = 0;

    while (status >= 0) {
        // The status at this level is just the error code.
        // It is recommended to define the error code as negative integer.
        status = ybim.evolve(logFile);
        
        // There is an error during the MD simulation
        if (status < 0) {
            cerr << "Error: main()\n"
            << "\tybim run-time error: evolve()" << status << endl;
            logFile.close();
            exit(1); // The main level error code.
        }

        // Let ybim calculates and save the post-evolution
        //  The variable status is used for counting the accumulated samples.
        // In-Gee Kim, December, 2015
        status = ybim.postEvolve(logFile);

        if (status < 0) {
            // There is an error during the post-evolution
            cerr << "Congratulations:\n"
                << "\tThe autocorrelation functions are converged."
                << status << endl;
            logFile.close();
            break; // The main level error code.
        }

        cout << "The sample id " << (status - 1)
            << " in the ensemble has been run."
            << endl << endl;

        if (status > N_samples) {
            cout << "\nThe ensemble average over autocorrelation functions"
                << " are not complete." << endl;
            break;
        }
    } // while (status >= 0)

    // All the jobs are done, close the logFile
    logFile.close();
    cout << endl << "ybim finished\n" << "#####" << endl;

    return 0;
}
