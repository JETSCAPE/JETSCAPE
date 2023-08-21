#ifndef THERMPTNSAMPLER_H
#define THERMPTNSAMPLER_H

#include "JetScapeLogger.h"
#include <vector>
#include <random>


using namespace Jetscape;

class ThermalPartonSampler
{
 public:
	//constructor - initializes variables and random number generator
	ThermalPartonSampler(unsigned int ran_seed);

	//sampler - samples thermal partons (u,d,s) from brick or 2+1d / 3+1d hydro
	void samplebrick();
	void sample_2p1d(double eta_max);
	void sample_3p1d(bool Cartesian_hydro=false);

	//to load a hypersurface
	void set_hypersurface(std::vector<std::vector<double>> surf_in){surface = surf_in;}

	//setters for params
	void brick_length_width(double len_bri, double wid_bri){L = 2.*len_bri + 4.; W = 2.*wid_bri + 4.; Time = len_bri;} // +4 gives 2fm additional brick in each direction
	void brick_flow(double vx_in, double vy_in, double vz_in){Vx = vx_in; Vy = vy_in; Vz = vz_in;}
	void brick_Tc(double brick_Tc){T = brick_Tc/GEVFM;}

	//getters for thermal partons
	int nTot(         ){return Plist.size();}
	int th_Evnum(int i){return (int)Plist[i][0];}
	int th_pid(  int i){return (int)Plist[i][1];}
	int th_orig( int i){return (int)Plist[i][2];}
	int th_stat( int i){return (int)Plist[i][11];}
	double th_px(int i){return      Plist[i][3];}
	double th_py(int i){return      Plist[i][4];}
	double th_pz(int i){return      Plist[i][5];}
	double th_e( int i){return      Plist[i][6];}
	double th_x( int i){return      Plist[i][7];}
	double th_y( int i){return      Plist[i][8];}
	double th_z( int i){return      Plist[i][9];}
	double th_t( int i){return      Plist[i][10];}
	int th_nL(){return num_ud;}
	int th_nS(){return num_s;}


 private:

	//thermal parton list
	std::vector<std::vector<double>> Plist;	// List of particles ( [0]->event number; [1]->particle ID; [2]->origin; [3-6]->Px,Py,Pz,E; [7-10]->x,y,z,t; [11]->Particle Status )
											// Same format as in shower data, event number is always 1, origin is always 0, particle status is always 0 (indicates a thermal quark)

	// Functions
	bool SplitSort (double goal, int floor, int ceiling, int quark);
	void MCSampler(double Temp, int quark);					//Samples momentum distribution, redefines NewP, NewX, NewY, NewZ variables - momentum in rest frame
	double FermiPDF (double Pc, double Mc, double Tc, double muc);		//Gives form of Fermi-Dirac distribution (no normalization factors!). Mu is currently given by 0 everywhere
	void CDFGenerator(double Temp, double M, int quark);				//generates cumulative distribution in thermal rest frame "quark" is 1 for light and 2 for strange

	// random number handling
	std::mt19937_64 rng_engine; //RNG - Mersenne Twist - 64 bit
	std::uniform_real_distribution<double> distribution{0.0, 1.0}; // Uniform distribution between 0 and 1
	// Function to generate a random number between 0 and 1
    double ran() {return distribution(rng_engine);}
	// Function to get the random number generator
    std::mt19937_64& getRandomGenerator() {return rng_engine;}

	// Static parameters, do not change
	const double PI = 3.141592653589793;
	double muPi0 = 0.;
	const double GEVFM = 0.197327053;

	// Adjustable parameters
	double xmq, xms, T, NUMSTEP;
	double CellDX, CellDY, CellDZ, CellDT;

	// Gaussian weights and abscissa (50 pt)
	// Used in numeric integration
	int GPoints = 50;		// Amount of Gaussian points for quadrature
	double GWeight[50];		// Gaussian weights for integration
	double GAbs[50];		// Gaussian abscissas for integration

	// Flags
	bool SetNum, SetNumLight, SetNumStrange, ShuffleList;

	// Global vars between functions
	double NewX, NewY, NewZ, NewP;
	int num_ud, num_s;
	std::vector<std::vector<double>> CDFTabLight;	// Cumulative distribution for thermal light quarks
	std::vector<std::vector<double>> CDFTabStrange;	// Cumulative distribution for s quarks

	// Brick info
	//L is the length of the brick sampled for thermal partons
	//W is the width of the face of the brick - should be scaled somewhat against L
	//Vx,Vy,Vz is a flow given to the brick
	double L, W, Time, Vx, Vy, Vz;

	// HyperSurface
	std::vector<std::vector<double>> surface;
};


#endif // THERMPTNSAMPLER_H