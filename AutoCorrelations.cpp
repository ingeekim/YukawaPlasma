//
// AutoCorrelations class implementation
//
// In-Gee Kim   December, 2015
//

// Private header
#include "AutoCorrelations.h"

AutoCorrelations::AutoCorrelations(int N_samples, int idx,
                                   unsigned long M_timeSteps,
                                   unsigned long M_snapShots,
                                   unsigned long N_particles)
    : M_steps(M_timeSteps / M_snapShots),
        N_types(2), N_ptls(N_particles),
        vxt(M_steps, vector<double>(N_ptls, 0.0)),
        vyt(M_steps, vector<double>(N_ptls, 0.0)),
        vzt(M_steps, vector<double>(N_ptls, 0.0)),
        jxt(M_steps, 0.0), jyt(M_steps, 0.0), jzt(M_steps, 0.0),
        Zv(N_types, vector<double>(M_steps, 0.0)),
        Zj(M_steps, 0.0),
        accZv(N_types, vector<double>(M_steps, 0.0)),
        accZj(M_steps, 0.0)
{
    //
    // AutoCorrelations class constructor
    //

    idxEnsemble = idx;

} // end of AutoCorrelations class constructor

AutoCorrelations::~AutoCorrelations()
{
    //
    // AutoCorrelations class destructor
    //

} // end of AutoCorrelations class destructor

void AutoCorrelations::setSpecies(unsigned long N_l, unsigned long N_h,
                                          double c_l, double c_h)
{
    N_light = N_l;
    N_heavy = N_h;
    x_l = c_l;
    x_h = c_h;
} // end of AutoCorrelations::set_individualParticles()

void AutoCorrelations::setUnits(double t_s, double r_s, double omega_p)
{
    half = (1.0 / 2.0);
    oneThird = (1.0 / 3.0);

    h = t_s;
    a_l = r_s;
    w_l = omega_p;
}
