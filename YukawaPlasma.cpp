//
// YukawaPlasma class implementation
//
// In-Gee Kim   July, 2015
//

// Private header
#include "YukawaPlasma.h"

YukawaPlasma::YukawaPlasma()
{
    //
    // YukawaPlasma class constructor
    //
} // end of YukawaPlasma class constructor

YukawaPlasma::YukawaPlasma(ofstream &logFile)
{
    //
    // YukawaPlasma class constructor
    //

    // input file open
    ifstream inputFile("input.yukawa", ios::in);
    if ( !inputFile ) {
        cerr << "Error: YukawaPlasma::YukawaPlasma\n"
            << "\tMissing input.yukawa file!" << endl;

        inputFile.close();
        logFile.close();
        exit(10); // error code at the YukawaPlasma constructor
    }

    // parse the input parameters
    int parseError = parseInput(inputFile);
    if ( parseError != 0 ) {
        cerr << "Error: YukawaPlasma::YukawaPlasma\n"
            << "\tParsing input.yukawa by parseInput() has problem(s)."
            << endl;
        inputFile.close();
        logFile.close();
        exit(10);   // error code at the YukawaPlasma constructor
    }

    inputFile.close();

    preparePlasma(logFile);

    cout << "A Yukawa binary ion mixture is created." << endl << endl;
} // end of YukawaPlasma class constructor

YukawaPlasma::~YukawaPlasma()
{
    //
    // YukawaPlasma class destructor
    //

} // end of YukawaPlasma class destructor

void YukawaPlasma::preparePlasma(ofstream &logFile)
{
    //
    // Set up the ensemble average index
    //
    
    idxEnsemble = 1;

    //
    // Prepare the plasma parameters from the input parameters
    //

    // set up the mathematical constants
    const double eps = 1.0e-12;
    pi = 4.0 * atan(1.0);
    twoPi = 2.0 * pi;
    fourPi = 4.0 * pi;
    inv_fourPi = 1.0 / fourPi;
    three = 3.0;
    half = 1.0 / 2.0;
    oneThird = 1.0 / 3.0;
    twoThird = 2.0 / 3.0;
    fiveThird = 5.0 / 3.0;
    threeHalves = 3.0 / 2.0;
    threePiSqtwoThird = pow((3.0 * pi * pi), twoThird);
    r_unitSphere = pow((3.0 / fourPi ), oneThird);

    // set up the physical constants
    E_H = 27.21138602;
    a_B = 0.52917721092e-8;
    amu = 931.494061e6;
    emu = 0.510998910e6;
    csq = 8.986153015714640000e20;

    esq = E_H * a_B;

    //
    // set up the plasma parameters
    //
    // The convention follows after
    // M. V. Beznogov and D. G. Yakovlev, Phys. Rev. E 90, 033102 (2014)
    //

    // ionic part
    N_particles = N_light + N_heavy;
    x_l = ((double) N_light) / ((double) N_particles);
    x_h = ((double) N_heavy) / ((double) N_particles);
    m_bar = x_l * m_light + x_h * m_heavy;

    a_l = r_unitSphere * pow(n_light, -oneThird); // reference length
    bar_a = a_l / a_B;
    if ( (N_heavy == 0) || (n_heavy < eps) ) {  // one-component simulation
        n_heavy = 0.0;
        a_h = 0.0;
    } else {
        a_h = r_unitSphere * pow(n_heavy, -oneThird);
    }

    n_ions = n_light + n_heavy;
    r_s = r_unitSphere * pow(n_ions, -oneThird);
    edgeL = pow((((double) N_particles) / n_ions), oneThird) / a_l;
    halfL = edgeL / 2.0;

    double T_ion = 0.0;
    if ( Gamma_light < 0.0 ) {
        // set up Gamma_0
        Gamma_0 = fabs(Gamma_light);
        T_ion = ((Z_light * Z_heavy * esq) / r_s) / Gamma_0;
        T_light = T_ion;
        // adjust Gamma_light
        Gamma_light_adjust = ((Z_light * Z_light * esq) / a_l) / T_light;
        // set Gamma_heavy = 0.0;
        Gamma_heavy = 0.0;
    } else {
        T_light = ((Z_light * Z_light * esq) / a_l) / Gamma_light;
    }

    if ( fabs(Gamma_heavy) < eps ) {
        T_heavy = T_light;
        Gamma_heavy_adjust = ((Z_heavy * Z_heavy * esq) / a_h) / T_heavy;
    } else {
        T_heavy = ((Z_heavy * Z_heavy * esq) / a_h) / Gamma_heavy;
    }

    if ( Gamma_light < 0.0 ) {
        Gamma_bar = x_l * Gamma_light_adjust + x_h * Gamma_heavy_adjust;
    } else {
        Gamma_bar = x_l * Gamma_light + x_h * Gamma_heavy;
    }

    // electronic part
    n_e = Z_light * n_light + Z_heavy * n_heavy;
    E_F = pow((3.0 * pi * pi * n_e), twoThird);
    E_F *= E_H * a_B * a_B * half;

    T_e = 0.0;
    if ( kappa < 0.0 ) {
        // the electronic temperature is set to T_ion
        T_e = T_ion;
        if ( quantumMode == 0 ) {
            // fully classical Debye length
            el_Debye = sqrt(T_e / (fourPi * n_e * esq));
        } else if (quantumMode == 1 ) {
            // Murillo  Debye length
            el_Debye = sqrt(
                        sqrt(T_e * T_e + (twoThird * E_F) * (twoThird * E_F))
                            / (fourPi * n_e * esq));
        } else {
            // no implementation
            cerr << "Error in input section &electron:" << endl;
            cerr << "\tThe quantumMode " << quantumMode << " is not supported.";
            cerr << endl;
            exit(11);
        }
        lambda_D = el_Debye / a_l;
        kappa_adjust = 1.0 / lambda_D;
        // kappa_adjust is assigned to kappa at the function
        // void YukawaPlasma::reportInputDiagnoses(ofstream &logFile)
    } else {
        // The normal case
        T_e = (fourPi * esq * n_e) * (a_l / kappa) * (a_l / kappa);
        lambda_D = 1.0 / kappa;
        el_Debye = lambda_D * a_l;

        // for quantum mode
        if ( quantumMode == 1 ) {
            // Murillo  Debye length
            el_Debye = sqrt(
                        sqrt(T_e * T_e + (twoThird * E_F) * (twoThird * E_F))
                            / (fourPi * n_e * esq));
            lambda_D = el_Debye / a_l;
            kappa_adjust = 1.0 / lambda_D;
            // kappa_adjust is assigned to kappa at the function
            // void YukawaPlasma::reportInputDiagnoses(ofstream &logFile)
        }
    }

    // equal temperature case
    Z_bar = n_e / n_ions;
    Z_barFiveThird = x_l * pow(Z_light, fiveThird);
    Z_barFiveThird += x_h * pow(Z_heavy, fiveThird);
    if ( Gamma_light > 0.0 ) {
        // When Gamma_light < 0.0, Gamma_0 is given by input
        Gamma_0 = Gamma_bar / (Z_barFiveThird * pow(Z_bar, oneThird));
    }

    // set up the plasma frequency
    w_l = sqrt((fourPi * n_light * Z_light * Z_light * esq * csq)
               / (m_light * amu));
    w_h = sqrt((fourPi * n_heavy * Z_heavy * Z_heavy * esq * csq)
               / (m_heavy * amu));

    w_e = sqrt((fourPi * n_e * esq * csq) / emu);
    w_Vlasov = sqrt(w_l * w_l + w_h * w_h);
    w_hydro = sqrt((fourPi * n_ions * Z_bar * Z_bar * esq * csq)
                   / (m_bar * amu));


    // Prepare the position, velocities, and acceleration arrays
    rx = new double[N_particles];
    ry = new double[N_particles];
    rz = new double[N_particles];
    vx = new double[N_particles];
    vy = new double[N_particles];
    vz = new double[N_particles];
    ax = new double[N_particles];
    ay = new double[N_particles];
    az = new double[N_particles];

    // Prepare the mass, charge, and temperature arrays
    mass = new double[N_particles];
    Z = new double[N_particles];
    T = new double[N_particles];

    Mass = 0.0;
    for (unsigned long i = 0; i < N_particles; i++) {
        if (i < N_light) {  // light components
            mass[i] = 1.0;
            Z[i] = 1.0;
            T[i] = T_light;
        } else { // heavy components
            mass[i] = (m_heavy / m_light);
            Z[i] = (Z_heavy / Z_light);
            T[i] = T_heavy;
        }
        Mass += mass[i];
    } // for (unsigned long i = 0; i < N_particles; i++)

    // Prepare the parameters for the pair correlation functions g_ij(r)
    R_bins = 500;
    r_max = 5.0;
    dr = r_max / ((double) R_bins);
    g_ll = new double[R_bins];
    g_hh = new double[R_bins];
    g_lh = new double[R_bins];
    for (unsigned long i = 0; i < R_bins; i++) {
        g_ll[i] = g_hh[i] = g_lh[i] = 0.0;
    } // for (unsigned long i = 0; i < R_bins; i++)

    reportInputDiagnoses(logFile);
}

void YukawaPlasma::reportInputDiagnoses(ofstream &logFile)
{
    const double eps = 1.0e-12;

    // write the information to the logFile
    logFile << "Number of light ions = " << N_light << endl;
    logFile << "Gamma of light ions = " << Gamma_light << endl;
    logFile << "Light ion density = " << n_light << " /cc" << endl;
    logFile << "Light ion temperature = " << T_light << " eV" << endl;
    logFile << "Light ion mass = " << m_light << " u" << endl;
    logFile << "Z_light = " << Z_light << endl;
    logFile << "Mole fraction of light ions = " << x_l << endl;
    logFile << "Plasma frequency of light ions = " << w_l << " Hz" << endl;
    logFile << endl;

    if ( Gamma_light < 0.0 ) {
        logFile << "The Gamma was set to total ionic Gamma_0: " << endl;
        logFile << "\tThe adjusted Gamma of light ions = "
            << Gamma_light_adjust << endl;
    }
    logFile << endl;

    logFile << "Number of heavy ions = " << N_heavy << endl;
    logFile << "Gamma of heavy ions = " << Gamma_heavy << endl;
    logFile << "Heavy ion density = " << n_heavy << " /cc" << endl;
    logFile << "Heavy ion temperature = " << T_heavy <<  " eV" << endl;
    logFile << "Heavy ion mass = " << m_heavy << " u" << endl;
    logFile << "Z_heavy = " << Z_heavy << endl;
    logFile << "Mole fraction of heavy ions = " << x_h << endl;
    logFile << "Plasma frequency of heavy ions = " << w_h << " Hz" << endl;
    logFile << endl;

    if ( fabs(Gamma_heavy) < eps ) {
        logFile << "The temperatures were set to the same:" << endl;
        logFile << "\tThe adjusted Gamma of heavy ions = "
            << Gamma_heavy_adjust << endl;
    }
    logFile << endl;


    // The total ionic information
    logFile << endl;
    logFile << "Total number of ions = " << N_particles << endl;
    logFile << "Total ion number density = " << n_ions << " /cc" << endl;
    logFile << endl;

    logFile << "Gamma_0 = " << Gamma_0 << endl;
    logFile << "The effective charge Z_bar = " << Z_bar << endl;
    logFile << "The number weighted effect charge Z_barFiveThird = "
            << Z_barFiveThird << endl;
    logFile << "Gamma_bar = " << Gamma_bar << endl;
    logFile << endl;
    logFile << "m_heavy / m_light = " << m_heavy / m_light << endl;
    logFile << "m_bar = " << m_bar << endl;
    logFile << "Z_heavy / Z_light = " << Z_heavy / Z_light << endl;
    logFile << endl;
    logFile << "The Vlasov ionic plasma frequency = " << w_Vlasov
        << " Hz" << endl;
    logFile << "                                  = " << w_Vlasov / w_l
        << " w_l" << endl;
    logFile << "The hydrodynamic ionic plasma frequency = " << w_hydro
        << " Hz" << endl;
    logFile << "                                        = " << w_hydro / w_l
        << " w_l" << endl;
    logFile << endl;

    logFile << "a_l = " << a_l << " cm" << endl;
    logFile << "a_h = " << a_h << " cm" << endl;
    logFile << "a_h / a_l = " << a_h / a_l << endl;
    logFile << "r_s = " << r_s << " cm" << endl;
    logFile << "r_s / a_l = " << r_s / a_l << endl;
    logFile << "Edge length in a_l = " << edgeL << endl;
    logFile << "half Edge length in a_l = " << halfL << endl;
    logFile << "Edge length in r_s = " << edgeL * (a_l / r_s) << endl;
    logFile << "half Edge length in r_s = " << halfL * (a_l / r_s) << endl;
    logFile << endl;

    // electronic part
    if ( (kappa < 0.0) || (quantumMode > 0) ) {
        logFile << "The electron temperature was adjusted" << endl;
        kappa = kappa_adjust;
        logFile << "\t kappa is adjusted." << endl;
    }
    logFile << "kappa = " << kappa << endl;
    logFile << "Electron density = " << n_e << " /cc" << endl;
    logFile << "Electron Fermi Energy = " << E_F << " eV" << endl;
    logFile << "Electron temperature = " << T_e << " eV" << endl;
    logFile << "quantumMode = " << quantumMode << endl;
    logFile << "Electron Debye length = " << lambda_D << " a_l" << endl;
    logFile << "                      = " << el_Debye << " cm" << endl;
    logFile << "Electronic plasma frequency = " << w_e << " Hz" << endl;
    logFile << endl;


    // operation part
    r2_cut = ((3.0 * 3.0)/ (kappa * kappa));
    logFile << "Radius square of the interaction cut-off = " << r2_cut << endl;
    logFile << "dt = " << dt << " /\\omega_l" << endl;
    logFile << "Number of steps = " << M_timeSteps << endl;
    logFile << "Steps per snapshots = " << M_snapShots << endl;
    logFile << "Number of equilibrium steps = " << M_eqlb << endl;
    logFile << "Number of pre-DIH steps = " << M_preDIH << endl;
    logFile << "Thermostat mode = " << thermostat << endl;

} // end of YukawaPlasma::reportInputDiagnoses()


void YukawaPlasma::set_ensemble(int Number_of_samples_in_the_Ensemble)
{
    N_samples = Number_of_samples_in_the_Ensemble;
}
