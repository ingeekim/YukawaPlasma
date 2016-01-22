//
// YukawaPlasma class implementation
// for analysis
//
// In-Gee Kim, August 2015
//

// Private header
#include "YukawaPlasma.h"


void YukawaPlasma::vel_dist(long mode)
{
    //
    // calculates the velocity distribution function
    //
    double v_max = 0.0;
    double *v = new double[N_particles];
    
    // finding the maximum velocity
    for (unsigned long i = 0; i < N_particles; i++) {
        v[i] = sqrt(vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
        if ( v[i] > v_max ) v_max = v[i];
    } // for (unsigned long i = 0; i < N_particles; i++)

    const unsigned long N_bins = 100;
    double dv = v_max / ((double) N_bins);

    double f[N_bins] = {0.0, };
    double f_l[N_bins] = {0.0, };
    double f_h[N_bins] = {0.0, };
    unsigned long igin = 0;
    for (unsigned long i = 0; i < N_particles; i++) {
        igin = (unsigned long) (v[i] / dv);
        if ( igin < N_bins ) {
            if ( i < N_light ) {
                f_l[igin] += 1.0;
            } else {
                f_h[igin] += 1.0;
            }
            f[igin] += 1.0;
        } // if (igin < N_bins)
    } // for (unsigned long i = 0; i < N_particles; i++)

    string fileName;
    if ( mode == 0 ) {
        // special mode for initial distribution
        fileName = "init_fv.out";
    } else if ( mode == -1 ) {
        // special mode for equilibrated distribution
        fileName = "eqlb_fv.out";
    }else if ( mode == -2 ) {
        // special mode for final distribution
        fileName = "final_fv.out";
    } else {
        // general mode for each loop
        fileName = "loop";
        string lp = to_string(mode);
        fileName += lp;
        fileName += "_fv.out";
    }

    ofstream distFile(fileName, ios::out);

    double vr;
    for (unsigned long i = 0; i < N_bins; i++) {
        vr = (((double) i) + half) * dv;
        f[i] /= N_particles;
        f_l[i] /= N_particles;
        f_h[i] /= N_particles;

        distFile << right << fixed << setw(15) << setprecision(9)
            << vr << " "
            << scientific << setw(24) << setprecision(15)
            << f[i] << " " << f_l[i] << " " << f_h[i] << endl;
    }

    delete v;
    distFile.close();

} // end of YukawaPlasma::init_Velocities()

