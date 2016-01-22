//
// AutoCorrelations Class definition
//
// In-Gee Kim, December 2015
//
#ifndef _AutoCorrelations
#define _AutoCorrelations

//
// System headers
//

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

#include <vector>
using std::vector;


// The AutoCorrelations class definition
class AutoCorrelations {
public:
    // the constructor with the member initializer
    AutoCorrelations(int, int, unsigned long, unsigned long, unsigned long);


    // destructors
    ~AutoCorrelations();

    // variable setting
    void setSpecies(unsigned long, unsigned long, double, double);
        // set N_light and N_heavy
    void setUnits(double, double, double); // set h, a_l, and w_l

    void readVelocities(void); // save the final condition
    void setCurrents(void); // calculate the current density
    void calcAutoCorr(void); // calculate the autocorrelation functions
    void accAutoCorr(void); // accumulates the autocorrelation functions
    void avgAutoCorr(void); // avergae the autocorrelation functions
    void savAutoCorr(void); // save the autocorrelation functions
    int chkStatFit(ofstream &); // check the statistical fitness
    void getDarken(ofstream &); // calculate the diffusion coefficients
                                //   by Darken Rule
    void getFick(ofstream &); // calculate the Fickian
                                // interdiffusion coefficient

private:
    //
    // private member constants
    //
    const unsigned long M_steps; // the number of time steps for autocorrelations
    const unsigned long N_types; // the number of types
    const unsigned long N_ptls;  // the number of particles

    // private member constant values
    unsigned long N_light; // the number of light particles
    unsigned long N_heavy; // the number of heavy particles
    double x_l; // the concentration of light species
    double x_h; // the concentration of heavy species
    double h; // the time step size
    double a_l; // the Wigner-Seitz radius of the light species
    double w_l; // the plasma frequency of the light species
    double half; // 1/2
    double oneThird; // 1/3

    // private member variables
    int N_samples;
    int idxEnsemble;

    //
    // private member vectors
    //
    vector<vector<double> > vxt;  // the x-component velocities
    vector<vector<double> > vyt;  // the y-component velocities
    vector<vector<double> > vzt;  // the z-component velocities

    vector<double> jxt;  // the x-component current
    vector<double> jyt;  // the y-component current
    vector<double> jzt;  // the z-component current

    vector<vector<double> > Zv;  // velocity autocorrelation functions
    vector<double> Zj; // current autocorrelation function

    vector<vector<double> > accZv; // accumulated vacf
    vector<double> accZj; // accumulated jacf

}; // end class YukawaPlasma

#endif
