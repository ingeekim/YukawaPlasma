//
// YukawaPlasma class implementation
//   spherical Yukawa force calculation
//
// In-Gee Kim   July, 2015
//

// Private header
#include "YukawaPlasma.h"

double YukawaPlasma::Ewald_parameter(void)
{
    // Calculates the necessary information of
    // the spherical approximation of the Ewald summation
    // for the Yukawa potential.
    // The approximation is appeared in Eq. (3.9)
    // J. M. Caillol and D. Gilles, J. Stat. Phys. Vol. 100, 905 (2000).
    //
    // In-Gee Kim, July 2015
    const int Max = 24;

    double length = 0.0;
    double dl, dm, dn;
    double alpha;

    double Ewald_factor = 0.0;
    for (int l = -Max; l <= Max; l++) {
        for (int m = -Max; m <= Max; m++) {
            for (int n = -Max; n <= Max; n++) {
                dl = (double) l;
                dm = (double) m;
                dn = (double) n;
                length = sqrt(dl * dl + dm * dm + dn * dn);

                if (length != 0.0) {
                    alpha = length * edgeL / lambda_D;
                    Ewald_factor += exp(-alpha) / (length * edgeL);
                } // if (length != 0.0)
            } // for (int n = -Max; n <= Max; n++)
        } // for (int m = -Max; m <= Max; m++)
    } // for (int l = -Max; l <= Max; l++)

    return Ewald_factor;
} // end of YukawaPlasma::Ewald_parameter()

double YukawaPlasma::spherical_Yukawa(void)
{
    //
    // Calculates the forces based on the yukawa potential.
    // Returns the potential energy.
    //
    // In-Gee Kim, July 2015

    const double halfL = edgeL / 2.0;  // half-length
    double potential_Energy = 0.0;

    // Get rid of old forces before summing.
    for (unsigned long i = 0; i < N_particles; i++) {
        ax[i] = ay[i] = az[i] = 0.0;
    }

    // define the inner loop variables for faster looping
    double dx, dy, dz = 0.0;
    double xj, yj, zj = 0.0;
    double dist, dist2, dist3 = 0.0;
    double kappa_r = 0.0;
    double U_ij, exp_mr, exp_pr, sinh_r, cosh_r = 0.0;
    double term1, term2, factor = 0.0;
    double m_r, dax, day, daz = 0.0;

    for (unsigned long i = 0; i < N_particles; i++) {
        for (unsigned long j = i + 1; j < N_particles; j++) {
            // obtain the primitive distance between two particles
            dx = rx[i] - rx[j];
            xj = rx[j];
            dy = ry[i] - ry[j];
            yj = ry[j];
            dz = rz[i] - rz[j];
            zj = rz[j];

            // minimum image shifts
            if ( dx < -halfL ) xj = rx[j] - edgeL;
            if ( dx > halfL ) xj = rx[j] + edgeL;
            if ( dy < -halfL ) yj = ry[j] - edgeL;
            if ( dy > halfL ) yj = ry[j] + edgeL;
            if ( dz < -halfL ) zj = rz[j] - edgeL;
            if ( dz > halfL ) zj = rz[j] + edgeL;

            // calculate the distance powers
            dx = rx[i] - xj;
            dy = ry[i] - yj;
            dz = rz[i] - zj;
            dist2 = dx * dx + dy * dy + dz * dz;
            if (dist2 > r2_cut) continue;

            dist = sqrt(dist2);
            dist3 = dist * dist2;

            //
            // Yukawa force and potential energy
            // under the spherical Ewald approximation
            // see YukawaPlasma::Ewald_parameter()
            //
            kappa_r = dist * kappa;
            exp_mr = exp(-kappa_r);
            exp_pr = exp(kappa_r);
            sinh_r = (exp_pr - exp_mr) / 2.0;
            cosh_r = (exp_pr + exp_mr) / 2.0;

            // Bare Yukawa term
            term1 = exp_mr * (1.0 + kappa_r);
            // Ewald summation term
            term2 = -(Ew / kappa) * (kappa_r * cosh_r - sinh_r);
            // total contribution
            factor = term1 + term2;
            factor *= (Z[i] * Z[j]) / (mass[i]);

            // add up to the potential energy
            U_ij = (exp_mr + (Ew / kappa) * sinh_r);
            U_ij *= ((Z[i] * Z[j] * Z_light * Z_light) / (dist * bar_a));
            U_ij *= E_H; // in eV
            potential_Energy += U_ij; // accumulate

            // force at current positions
            dax = (dx / dist3) * factor * oneThird;
            day = (dy / dist3) * factor * oneThird;
            daz = (dz / dist3) * factor * oneThird;

            // update the force
            m_r = mass[i] / mass[j];
            ax[i] += dax;
            ay[i] += day;
            az[i] += daz;
            ax[j] -= (dax * m_r);
            ay[j] -= (day * m_r);
            az[j] -= (daz * m_r);

        } // for (int j = i + 1; j < N_particles; j++)
    } // for (int i = 0; i < N_particles; i++)

    return potential_Energy;
} // end of YukawaPlasma::spherical_Yukawa()
