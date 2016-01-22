//
// YukawaPlasma class implementation
//    Thermostats
//
// In-Gee Kim   August, 2015
//

// Private header
#include "YukawaPlasma.h"

void YukawaPlasma::strong_scale(vector<double> T_current)
{
    double factor = 0.0;
    double T_cur = 0.0;

    for (unsigned long i = 0; i < N_particles; i++) {
        if ( i < N_light ) {
            T_cur = T_current[0];
        } else {
            T_cur = T_current[1];
        }
        factor = sqrt(T[i] / T_cur);
        vx[i] *= factor;
        vy[i] *= factor;
        vz[i] *= factor;
    }
} // end of YukawaPlasma::sclae_vel1()

void YukawaPlasma::weak_scale(vector<double> T_current)
{
    double factor = 0.0;
    double T_cur = 0.0;

    for (unsigned long i = 0; i < N_particles; i++) {
        if ( i < N_light ) {
            T_cur = T_current[0];
        } else {
            T_cur = T_current[1];
        }
        factor = sqrt((20.0 + T[i] / T_cur - 1.0) / 20.0);
        vx[i] *= factor;
        vy[i] *= factor;
        vz[i] *= factor;
    }
} // end of YukawaPlasma::sclae_vel1()

void YukawaPlasma::light_scale(vector<double> T_current)
{
    // This is a thermostat mode for test case.
    double factor = 0.0;
    double T_cur = 0.0;

    for (unsigned long i = 0; i < N_particles; i++) {
        if ( i < N_light) {
            T_cur = T_current[0];
            factor = sqrt(T[i] / T_cur);
            vx[i] *= factor;
            vy[i] *= factor;
            vz[i] *= factor;
        }
    }
} // end of YukawaPlasma:::test_scale()

void YukawaPlasma::heavy_scale(vector<double> T_current)
{
    // This is a thermostat mode for test case.
    double factor = 0.0;
    double T_cur = 0.0;

    for (unsigned long i = 0; i < N_particles; i++) {
        if ( i > N_light) {
            T_cur = T_current[1];
            factor = sqrt(T[i] / T_cur);
            vx[i] *= factor;
            vy[i] *= factor;
            vz[i] *= factor;
        }
    }
} // end of YukawaPlasma:::test_scale()
