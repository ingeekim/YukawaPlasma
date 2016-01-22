//
// YukawaPlasma class implementation
//
// In-Gee Kim   July, 2015
//

// Private header
#include "YukawaPlasma.h"

vector<double> YukawaPlasma::velocity_Verlet(double &potential_Energy)
{
    //
    // Advance the particle positions and velocities.
    // The particle posisions and velocities are computed
    // at the next time step using teh velocity Verlet algorithm.
    //
    // The argument potential_Energy is the potential energy returns.
    //
    // Returns the kinetic energy.
    //
    // In-Gee Kim, July 2015

    // define the loop variables outside the loop for speed-up
    double dx, dy, dz =0.0;
    const double halfL = edgeL / 2.0;

    // first phase advancement
    for (unsigned long i = 0; i < N_particles; i++) {
        // Compute amounts of position advancements
        dx = vx[i] * dt + half * ax[i] * dt * dt;
        dy = vy[i] * dt + half * ay[i] * dt * dt;
        dz = vz[i] * dt + half * az[i] * dt * dt;

        // Advance the positions
        rx[i] += dx;
        ry[i] += dy;
        rz[i] += dz;

        // Regulate the minimum image shifts
        if ( rx[i] < -halfL ) rx[i] += edgeL;
        if ( rx[i] > halfL ) rx[i] -= edgeL;
        if ( ry[i] < -halfL ) ry[i] += edgeL;
        if ( ry[i] > halfL ) ry[i] -= edgeL;
        if ( rz[i] < -halfL ) rz[i] += edgeL;
        if ( rz[i] > halfL ) rz[i] -= edgeL;

        // Advance the velocities
        vx[i] += half * ax[i] * dt;
        vy[i] += half * ay[i] * dt;
        vz[i] += half * az[i] * dt;
    } // for (int i = 0; i < N_particles; i++)

    // Obtain new forces
    potential_Energy = spherical_Yukawa();

    // Update velocities to t + dt
    //    using intermediate velocities and new forces.
    vector<double> kinetic_Energy(2, 0.0);

    double kinetic_i = 0.0;
    for (unsigned long i = 0; i < N_particles; i++) {
        vx[i] += half * ax[i] * dt;
        vy[i] += half * ay[i] * dt;
        vz[i] += half * az[i] * dt;
        kinetic_i = vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
        kinetic_i *= threeHalves * mass[i] * Z_light * Z_light;
        kinetic_i *= (E_H / bar_a); // in eV
        if ( i < N_light ) {
            kinetic_Energy[0] += kinetic_i; // accumulate light component
        } else {
            kinetic_Energy[1] += kinetic_i; // accumulate heavy component
        } // if ( i < N_light )
    } // for (unsigned long i = 0; i < N_particles; i++)

    return kinetic_Energy;
} // end of YukawaPlasma::velocity_Verlet()
