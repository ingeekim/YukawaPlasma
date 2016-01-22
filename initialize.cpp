//
// YukawaPlasma class implementation
// for the Plasma initialization
//
// In-Gee Kim, July 2015
//

// Private header
#include "YukawaPlasma.h"

void YukawaPlasma::init_Positions(void)
{
    // Initialize Positions.
    // This is based on the Fortran routine of yukawa.f90 by M. S. Murillo.
    // The particle positions are initialized randomly in a cube
    // with lengths measured in ion-sphere radii units (r_s).
    // Particles are not placed too close to the edge of the cube
    // nor they placed too close together.
    // The cloesest two particles may be is one ion-sphere radius.
    //
    // In-Gee Kim, July 2015

    // random number generation
    srand(time(0));

    const unsigned long NT = 1000000;
    const double dNT = (double) NT;

    // Assign a position to each particle.
    double dist, xx, yy, zz = 0.0;
    double shortest = 0.01;
    bool is_Overlap;
    for(unsigned long i = 0; i < N_particles; i++) {
        is_Overlap = true;
        // While particles overlap aonther, keep choosing new random positions.
        unsigned long count_overlap = 0;
        while (is_Overlap) {
            // Assume it doesn't overlap with other particles already assigned.
            is_Overlap = false;
            // Pick a position not within 0.5 ion-sphere radii
            // of the cube's edge.
            xx = ((double) (rand() % NT) / dNT) * (edgeL - 1.0) + 0.5;
            yy = ((double) (rand() % NT) / dNT) * (edgeL - 1.0) + 0.5;
            zz = ((double) (rand() % NT) / dNT) * (edgeL - 1.0) + 0.5;
            // Check to see if this is within one ion-sphere radius
            // of any particle.
            for(unsigned long j = 0; j < i; j++) {
                dist = (rx[j] - xx) * (rx[j] - xx)
                + (ry[j] - yy) * (ry[j] - yy) + (rz[j] - zz) * (rz[j] - zz);
                if ( dist <= shortest ) {
                    is_Overlap = true;
                } // if-else
            } // for(j = 0; j < i; j++)
            count_overlap++;
            if (count_overlap > 100) {
                cerr << "The minimum particle distance is too long." << endl;
                exit(20);
            }
        } // while (is_Overlap)
        // Good position found. Store position.
        rx[i] = xx - halfL;
        ry[i] = yy - halfL;
        rz[i] = zz - halfL;
    } // for(i = 0; i < N_particles; i++)
} // end of YukawaPlasma::init_Positions()

void YukawaPlasma::init_Velocities(void)
{
    // Initialize Velocities.
    // This is based ont eh Fortran routine of yukawa.f90 by M. S. Murillo.
    // The particle velocities are initialized randomly
    // using the Box-Muller method.
    // This procedure yields a Gaussian velocity for each velocity component
    // of each particles.
    // Velocities are then shifted to ensure that the total linear momentum
    // is zero for possible use as a diagnostic.
    //
    // In-Gee Kim, July 2015

    // random number generation
    srand(time(0));

    // Scale velocity distribution to ion temperature in eV
    // with units of r_s * \omega_p
    // The Maxwell-Boltzmann distribution function
    // is exp(- \frac{mv^2}{2 k_B T}).
    // The argument of the exponential function should be scaled as
    // \frac{mv^2}{2 k_B T} * (\frac{r_s \omega_p}{r_s \omega_p})^2
    // The electron charge square e^2 is defined as
    // the Hartree times Bohr radius in units of Ha[ev] a_B [cm]
    // cf. Ha := e^2 / a_B
    //
    // In-Gee Kim, July 2015

    const unsigned long NT = 1000000;
    const double dNT = (double) NT;

    double sigma = 0.0;

    // variables for checking the total linear momentum to be zero
    double sum_px = 0.0;
    double sum_py = 0.0;
    double sum_pz = 0.0;

    // assign the velocities in a random manner
    // with the scaling of Maxwell-Boltzmann distribution function
    double r1, r2 = 0.0;
    double Gamma1i = 0.0;
    for (unsigned long i = 0; i < N_particles; i++) {
        Gamma1i = ((Z_light * Z_light * esq) / a_l) / T[i];
        sigma = 1.0 / sqrt(3.0 * mass[i] * Gamma1i);
        r1 = ((double) (rand() % NT) / dNT);
        r2 = ((double) (rand() % NT) / dNT);
        vx[i] = sigma * sqrt(-2.0 * log(r1)) * cos(twoPi * r2);
        r1 = ((double) (rand() % NT) / dNT);
        r2 = ((double) (rand() % NT) / dNT);
        vy[i] = sigma * sqrt(-2.0 * log(r1)) * cos(twoPi * r2);
        r1 = ((double) (rand() % NT) / dNT);
        r2 = ((double) (rand() % NT) / dNT);
        vz[i] = sigma * sqrt(-2.0 * log(r1)) * cos(twoPi * r2);

        // sum the momenta
        sum_px += mass[i] * vx[i];
        sum_py += mass[i] * vy[i];
        sum_pz += mass[i] * vz[i];
    } // for(int i = 0; i < N_paricles; i++)

    // enforce the total linear momentum to be zero
    for (unsigned long i = 0; i < N_particles; i++) {
        vx[i] -= sum_px / Mass;
        vy[i] -= sum_py / Mass;
        vz[i] -= sum_pz / Mass;
    } // for(int i = 0; i < N_paricles; i++)

    // calculates and save the initial velocity distribution
    // In-Gee Kim, August 2015
    vel_dist(0);
} // end of YukawaPlasma::init_Velocities()

