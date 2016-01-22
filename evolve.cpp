//
// YukawaPlasma class implementation
//   evolution control
//
// In-Gee Kim   July, 2015
//

// Private header
#include "YukawaPlasma.h"

int YukawaPlasma::evolve(ofstream &logFile)
{
    //
    // YukawaPlasma class member function evolve()
    //
    // Invokes the molecular dynamics simulation indeed.
    //
    // It returns an int value 0 if the running has been successfully done,
    //   otherwise a certain error code, which will be defined later.
    //

    Ew = Ewald_parameter();
    logFile << endl;
    logFile << "Madelung constant --> " << Ew << endl;
    logFile << "Ewald factor --> " << (Ew / kappa) << endl;

    double potential_Energy, total_Energy = 0.0;
    vector<double> kinetic_Energy(2, 0.0);
    vector<double> T_current(2, 0.0);

    if ( initMode == 0 ) {
        logFile << "\nGenerate the system randomly" << endl;
        // generate the initial position and vectors
        init_Positions();
        logFile << "\tpositions are generated." << endl;
        init_Velocities();
        logFile << "\tvelocities are generated." << endl;

        // find forces at t = 0 to initialize acceleration array
        potential_Energy = spherical_Yukawa();

        // pre-equilibration Disorder-induced heating phase
        cout << "Pre-equilibration Disorder-induced heating ... "
            << M_preDIH << " steps" << endl;
        for (unsigned long i = 0; i < M_preDIH; i++) {
            kinetic_Energy = velocity_Verlet(potential_Energy);
        }

        // Equlibration phase
        cout << "Equilibrating ... " << M_eqlb << " steps" << endl;
        for (unsigned long i = 0; i < M_eqlb; i++) {
            kinetic_Energy = velocity_Verlet(potential_Energy);
            T_current[0] = kinetic_Energy[0] / ((double) N_light);
            T_current[0] *= twoThird;
            if ( N_heavy > 0 ) {
                T_current[1] = kinetic_Energy[1] / ((double) N_heavy);
                T_current[1] *= twoThird;
            } // if ( N_heavy > 0 )

            // scale velocites along with thermostat functions
            if ( thermostat == -4 ) {
                heavy_scale(T_current);
            } else if ( thermostat == -3 ) {
                light_scale(T_current);
            } else if ( thermostat == -2 ) {
                weak_scale(T_current);
            } else if ( thermostat == -1 ) {
                strong_scale(T_current);
            } else if ( thermostat != 0 ) {
                cout << "No thermostat mode " << thermostat << endl;
                exit(-1);
            }
        } // for (unsigned long i = 0; i < M_eqlb; i++)
    } else {
        // read the initial position and vectors
        cout << "Reading in initial conditions..." << endl;
        ifstream initialFile("initca.out", ios::in);

        for (unsigned long i = 0; i < N_particles; i++) {
            initialFile >> vx[i] >> vy[i] >> vz[i] >> rx[i] >> ry[i] >> rz[i];
            if ( initMode == 1 ) {  // accepting the old format
                rx[i] -= halfL;
                ry[i] -= halfL;
                rz[i] -= halfL;
            }
        } // for (unsigned long i = 0; i < N_particles; i++)

        initialFile.close();
    } // if ( initMode == 0 ) ... else ...

    // report the equilibriated velocity distributions
    // In-Gee Kim, August 2015
    vel_dist(-1);

    //
    // Main loop
    //

    // open files
    ofstream energyFile("evst.out", ios::out);
    ofstream temperatureFile("temp.out", ios::out);
    ofstream qpFile("iqp.out", ios::out);
    ofstream posFile("pos.xyz", ios::out);

    cout << "Entering main loop... " << M_timeSteps << " steps" << endl;
    for (unsigned long loop = 0; loop < M_timeSteps; loop++) {
        // advance particles
        kinetic_Energy = velocity_Verlet(potential_Energy);

        // calculate temperature
        T_current[0] = kinetic_Energy[0] / ((double) N_light);
        T_current[0] *= twoThird;
        if ( N_heavy > 0 ) {
            T_current[1] = kinetic_Energy[1] / ((double) N_heavy);
            T_current[1] *= twoThird;
        }

        // apply thermostats
        if ( thermostat == 1 ) {
            strong_scale(T_current);
        } else if ( thermostat == 2 ) {
            weak_scale(T_current);
        } else if ( thermostat == 3 ) {
            light_scale(T_current);
        } else if ( thermostat == 4 ) {
            heavy_scale(T_current);
        } else if ( thermostat > 4 ) {
            cout << "No thermostat mode for " << thermostat << endl;
            exit(-1);
        }

        // Every snapshot write outputs, find temp
        if ( loop % M_snapShots == 0 ) {
            //
            // writing snapt shot information
            //

            // temperature information
            temperatureFile << right << fixed << setw(15) << setprecision(9)
                << ((double) loop) * dt << " "
                << scientific << setw(24) << setprecision(15)
                << T_current[0] << " " << T_current[1] << endl;

            // energy information
            total_Energy = kinetic_Energy[0] + potential_Energy;
            if ( N_heavy > 0 ) {
                total_Energy += kinetic_Energy[1];
            }

            energyFile << right << fixed << setw(15) << setprecision(9)
                << ((double) loop) * dt << " "
                << scientific << setw(24) << setprecision(15)
                << total_Energy << " "
                << kinetic_Energy[0] << " " << kinetic_Energy[1] << " "
                << potential_Energy << endl;

            // velocities, positions, and accelerations
            posFile << right << fixed << setw(15) << N_particles << endl;
            posFile << "Lattice = " << "\""
                << fixed << setw(15) << setprecision(9)
                << edgeL << " " << "0.0  0.0 "
                << "0.0 " << edgeL << "0.0 "
                << "0.0 0.0 " << edgeL << "\"" << endl;
            for (unsigned long i = 0; i < N_particles; i++) {
                qpFile << right << scientific << setw(24) << setprecision(15)
                    << vx[i] << " " << vy[i] << " " << vz[i] << " "
                    << rx[i] << " " << ry[i] << " " << rz[i] << " "
                    << ax[i] << " " << ay[i] << " " << az[i] << endl;
                posFile << right << fixed << setw(15) << setprecision(9)
                    << Z[i] << " " << mass[i] * m_light << " "
                    << scientific << setw(24) << setprecision(15)
                    << rx[i] << " " << ry[i] << " " << rz[i] << " " << endl;
            } // for (unsigned long i = 0; i < N_particles; i++)

            // calculate the raidal distribution functions on the fly
            //  the g_ij(r) should be normalized before they are saved
            //
            // In-Gee Kim, August 2015
            double r, dx, dy, dz = 0.0;
            unsigned long igin = 0;
            for (unsigned long i = 0; i < N_particles; i++) {
                for (unsigned long j = 0; j < N_particles; j++) {
                    if ( i != j ) {
                        // distance
                        dx = rx[i] - rx[j];
                        dy = ry[i] - ry[j];
                        dz = rz[i] - rz[j];

                        // distance regularization
                        if ( dx > halfL ) dx -= edgeL;
                        if ( dx < -halfL ) dx += edgeL;
                        if ( dy > halfL ) dy -= edgeL;
                        if ( dy < -halfL ) dy += edgeL;
                        if ( dz > halfL ) dz -= edgeL;
                        if ( dz < -halfL ) dz += edgeL;

                        // radial distance
                        r = sqrt(dx * dx + dy * dy + dz * dz);
                        igin = (unsigned long) (r / dr);
                        if ( igin < R_bins ) {
                            if ( (i < N_light) && (j < N_light) ) {
                                g_ll[igin] += 1.0;
                            } else if ( (i > N_light) && (j > N_light) ) {
                                g_hh[igin] += 1.0;
                            } else {
                                g_lh[igin] += 1.0;
                            }
                        } // if ( igin < R_bins )
                    } // if ( i != j )
                } // for (unsinged long j = 0; j < N_particles; j++)
            } // for (unsigned long i = 0; i < N_particles; i++)
            // block end radial distribution function

        } // if ( loop % M_snapShots == 0 )

    } // for (unsigned long loop = 0; loop < M_timeSteps; loop++)

    // close files
    energyFile.close();
    temperatureFile.close();
    qpFile.close();
    posFile.close();

    return 0;
} // end of YukawaPlasma::evolve()
