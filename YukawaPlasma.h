//
// YukawaPlasma Class definition
//
// In-Gee Kim, July 2015
//
#ifndef _YukawaPlasma
#define _YukawaPlasma

//
// System headers
//
#include <cmath>
using std::atan;
using std::cos;
using std::exp;
using std::fabs;
using std::log;
using std::pow;
using std::sqrt;

#include <ctime>
using std::time;

#include <iomanip>
using std::fixed;
using std::right;
using std::scientific;
using std::setprecision;
using std::setw;

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using std::ios;

#include <fstream>
using std::ifstream;
using std::ofstream;

#include <cstdlib>
using std::exit;
using std::rand;
using std::srand;

#include <queue>
using std::queue;

#include <string>
using std::string;
using std::stoul;
using std::stod;
using std::to_string;

#include <vector>
using std::vector;

// The YukawaPlasma class definition
class YukawaPlasma {
public:
    // constructors
    YukawaPlasma();  // default constructor does nothing
    YukawaPlasma(ofstream &);  // this constructor prepares the system

    // destructors
    ~YukawaPlasma();

    int evolve(ofstream &);  // The instance of YukawaPlasma evolve in time
    int postEvolve(ofstream &); // The post-evolution processor

    double Ewald_parameter(void);  // The Ewald parameter calculator
    void init_Positions(void);  // Initialize the positions
    void init_Velocities(void); // Initialize the velocities
    double spherical_Yukawa(void); // calculate the force
    vector<double> velocity_Verlet(double &); // advacne particles
    void strong_scale(vector<double>); // strong velocity scale
    void weak_scale(vector<double>); // weak velocity scale
    void light_scale(vector<double>); // test velocity scale light only
    void heavy_scale(vector<double>); // test velocity scale heavy only

    void savFinal(ofstream &); // save the final condition
    void vel_dist(long); // velocity distribution
    void rad_dist(ofstream &); // save the radial distribution function
    void calcDiffusion(ofstream &); // calculates the diffusion coefficients

    void set_ensemble(int); // set upt the number of samples in the ensemble

private:
    //
    // private member function(s)
    //

    // the parser functions
    int parseInput(ifstream &);  // the input file parser
    int get_parseMode(string &); // determines the parser mode
    void parseLight(queue<string> &); // the light component parser
    void parseHeavy(queue<string> &); // the heavy component parser
    void parseEvolution(queue<string> &); // the evolution parser
    void parseElectron(queue<string> &); // the electron parser
    void parseOptions(queue<string> &); // the option parser
    void preparePlasma(ofstream &); // prepare the plasma parameters
    void reportInputDiagnoses(ofstream &); // reports input parse results

    //
    // private member variable(s)
    //

    // main function argument variable(s)
    int N_samples;  // the number of samples in ensemble

    // &heavy and &light sections
    unsigned long N_particles;   // the number of particles
    unsigned long N_light; // the light component
    unsigned long N_heavy; // the heavy component
    double Gamma_light;  // the Gamma parameter for the light
    double Gamma_light_adjust; // the light ion Gamma paramter adjusted
    double Gamma_heavy; // the Gamma parameter for the heavy
    double Gamma_heavy_adjust; // the heavy ion Gamma parameter adjusted
                         // when T_heavy = T_light
    double Gamma_bar; // the effective Gamma parameter
    double Gamma_0; // = Gamma_bar / (Z_bar^{5/3) Z_bar^{1/3})
    double n_light; // the light component number density in /cc
    double n_heavy; // the heavy component number density in /cc
    double n_ions; // the ion component number density in /cc
    double x_l; // the mole fraction of light component
    double x_h; // the mole fraction of heavy component
    double Z_light; // the charge of light component
    double Z_heavy; // the charge of heavy component
    double Z_bar; // the effective charge of the system
    double Z_barFiveThird; // the number weighted charge of the system
    double m_light; // the mass of light component
    double m_heavy; // the mass of heavy component
    double Mass; // the total normalized mass = \sum_i mass[i] / m_light

    // &evolution section
    double dt; // the time step size
    unsigned long M_timeSteps; // the number of time steps
    unsigned long M_snapShots; // the number of snap shots
    unsigned long M_eqlb;  // the number of equilibration time steps
    unsigned long M_preDIH; // the number of steps for DIH treatment before
                            // equilibration

    // &electron section
    double kappa; // the electron screening parameter
    double kappa_adjust; // the adjusted electron screening parameter
    double lambda_D;  // the electron Debye length in r_s
    double el_Debye;  // the electron Debye length in cm
    double T_e;  // the electron temperature in eV
    unsigned long quantumMode; // the mode for quantum effects of electrons
    double E_F; // the Fermi energy

    // &option section
    unsigned long initMode;  // the initial state mode
    unsigned long finalMode; // the final state mode
    long thermostat; // the main routine thermostat mode
    int idxEnsemble; // the index in the ensemble

    // the mathematical constants
    double pi;  // \pi
    double twoPi; // 2\pi
    double fourPi;  // 4\pi
    double inv_fourPi; // 1 / 4\pi
    double three; // 3
    double half; // 1/2
    double oneThird; // 1/3
    double twoThird; // 2/3
    double fiveThird; // 5/3
    double threeHalves; // 3/2
    double threePiSqtwoThird; // (3 \pi^2)^{2/3}
    double r_unitSphere; // (3 / 4\pi)^{1/3}

    // the physical constants
    double E_H;  // the Hartree energy in eV
    double a_B; // the Bohr radius in cm
    double esq; // the electron charge square e^2 in E_H * a_B
    double csq; // the square of speed of light in cm^2 / s^2
    double amu; // the atomic mass in eV / c^2
    double emu; // the electron mass in eV / c^2

    // the plasma parameters
    double n_e;  // the electron number density in /cc
    double m_bar; // the concentration weighted ionic mass
    double edgeL; // the edge length of the cell in cm
    double halfL; // the half edge length
    double T_light;  // the light ion component temperature in eV
    double T_heavy; // the heavy ion component temperature in eV
    double a_l;  // the Wigner-Seitz radius of the light component in cm
    double a_h;  // the Wigner-Seitz radius of the heavy component in cm
    double bar_a; // a_l / a_B
    double r_s;  // the Wigner-Seitz radius of the ion components in cm
    double Ew; // the Ewald parameter of hyperspherical approximation
    double w_l; // the light component plasma frequency in Hz
    double w_h; // the heavy component plasma frequency in Hz
    double w_e; // the electronic plasma frequency in Hz
    double w_Vlasov; // the ionic Vlasov plasma frequency in Hz
    double w_hydro; // the ionic hydrodynamic plasma frequency in Hz

    // Positions, velocities, accelerations, mass, and charges
    // units are normalized by the Wigner-Seitz radius a_1 for length
    // and by the plasma frequency 1/\omega_1 for time
    // of those of the light component.
    // The mass are normalized by the light component values.
    double *rx;
    double *ry;
    double *rz;
    double *vx;
    double *vy;
    double *vz;
    double *ax;
    double *ay;
    double *az;
    double *mass;  // mass ratio
    double *Z;  // charge
    double *T;  // temperature

    // Yukawa cutoff parameters
    double r2_cut;  // radius square of the interaction cut-off 

    // pair correlation function parameters
    unsigned long R_bins;
    double r_max;
    double dr;
    double *g_ll;
    double *g_lh;
    double *g_hh;

}; // end class YukawaPlasma

#endif
