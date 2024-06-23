#ifndef THERMPTNSAMPLER_H
#define THERMPTNSAMPLER_H

#include "JetScapeLogger.h"
#include <vector>
#include <random>
#include <omp.h>

using namespace Jetscape;

class ThermalPartonSampler
{
 public:
	//constructor - initializes variables and random number generator
	ThermalPartonSampler(unsigned int ran_seed, double hydro_Tc);

	//sampler - samples thermal partons (u,d,s) from brick or 2+1d / 3+1d hydro
	void samplebrick();
	void sample_2p1d(double eta_max);
	void sample_3p1d(bool Cartesian_hydro=false);

	//to load a hypersurface
	void set_hypersurface(std::vector<std::vector<double>> surf_in){surface = surf_in;}

	//setters for params
	void brick_length_width(double len_bri, double wid_bri){L = 2.*len_bri + 4.; W = 2.*wid_bri + 4.; Time = len_bri;} // +4 gives 2fm additional brick in each direction
	void brick_flow(double vx_in, double vy_in, double vz_in){Vx = vx_in; Vy = vy_in; Vz = vz_in;}

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

	// List of thermal partons sampled
	std::vector<std::vector<double>> Plist;
	// List of particles ( [0]->event number; [1]->particle ID; [2]->origin; 
	// [3-6]->Px,Py,Pz,E; [7-10]->x,y,z,t; [11]->Particle Status )
	// Same format as in shower data, event number is always 1, 
	// origin is always 0, particle status is always 0 (indicates a thermal quark)

	// random number handling
	std::mt19937_64 rng_engine; //RNG - Mersenne Twist - 64 bit
	std::uniform_real_distribution<double> distribution{0.0, 1.0}; // Uniform distribution between 0 and 1
	// Function to generate a random number between 0 and 1
    double ran() {return distribution(rng_engine);}
	// Function to get the random number generator
    std::mt19937_64& getRandomGenerator() {return rng_engine;}

	
	// HyperSurface
	std::vector<std::vector<double>> surface;

	// Vector to pre-tabulate the Bessel function K2
	std::vector<double> BesselK2;
	void InitializeBesselK2();
	// Bessel function K2
	double BesselK2function(double arg);

	// Fermi-Dirac distribution momentum integral
	// This uses the same method as iSS to calculate the momentum integral
	double FermiDiracDistributionMomentumIntegral(double T, double m);
	
	// Fermi-Dirac distribution without normalization factors
	double FermiDiracDistribution(double p, double m, double T);

	// Cumulative Distribution Function (CDF) generator in thermal rest frame  
	// (quark=1: light, quark=2 strange)
	void CDFGenerator(double T, double m, int quark);
	
	// Samples parton momentum in rest frame
	// Input temperature is assumed to be in fm^-1
	std::tuple<double, double, double, double> MomentumSampler(double T, int quark);

	// Set up the Lorentz boost matrix and the (gamma, v) vector
    void LorentzBoostMatrix(std::vector<double>& v, std::vector<std::vector<double>>& BoostMatrix, bool brick);

    // Perform Lorentz boost on vector
    std::vector<double> LorentzBoost(std::vector<double>& x, std::vector<std::vector<double>>& BoostMatrix);

    // Sample partons and fill Plist
    void SamplePartons(int Npartons, int quark, double T, bool brick,
		std::vector<double>& CPos, std::vector<std::vector<double>>& BoostMatrix,
		bool slice_boost, double eta_slice, std::vector<std::vector<double>>& Plocal);

	// Function to get the closest cached temperature to the target temperature - CDF tabulated
    double getClosestCachedTemp(const std::unordered_map<double, std::vector<std::vector<double>>>& cache, double targetTemp) const;
	// Function to create the CDF cache within a certain temperature range of five percent
	// using 20 temperature points in that range
	void createCDFCacheEntry(double& TRead, double xmq, int quarkType, 
                        std::unordered_map<double,
						std::vector<std::vector<double>>>& cache, 
                        std::vector<std::vector<double>>& CDFTab);
	// Function to get the closest cached temperature to the target temperature - Fermi-Dirac integral tabulated
	double getClosestCachedTempFermiDiracIntegral(const std::unordered_map<double, double>& cache, double targetTemp) const;

	// Constants
	double degeneracy_ud; // Degeneracy of UD quarks
	double degeneracy_s; // Degeneracy of S quarks
	double CDFcacheAccuracy;	// Accuracy range for cached temperature
	double xmq; // use pythia values here
	double xms; // use pythia values here
	int NUMSTEP; // 2^12+1, for steps of CDF Table, changes coarseness of momentum sampling

	int num_ud, num_s;	// Number of light and strange quarks sampled

	// Adjustable params for 3+1d
	// these are the values used in SurfaceFinder.cc
	double CellDX;
	double CellDY;
	double CellDZ;
	double CellDT;

	// Flags
	bool SetNum, SetNumLight, SetNumStrange, ShuffleList;

	std::vector<std::vector<double>> CDFTabLight; // Cumulative distribution for thermal light quarks
	std::vector<std::vector<double>> CDFTabStrange; // Cumulative distribution for s quarks

	// precompute CDFs and store them in a cache
	std::unordered_map<double, std::vector<std::vector<double>>> CacheCDFLight;
	std::unordered_map<double, std::vector<std::vector<double>>> CacheCDFStrange;

	// precompute part of the integral over the Fermi-Dirac distribution
	std::unordered_map<double, double> CacheFermiDiracIntegralLight;
	std::unordered_map<double, double> CacheFermiDiracIntegralStrange;

	// Brick parameters
	// L is the length of the brick sampled for thermal partons
	// W is the width of the face of the brick - should be scaled somewhat against L
	// Vx,Vy,Vz is a flow given to the brick
	double L, W, Time, Vx, Vy, Vz;
	double hydroTc; // converted from GeV into fm^-1
};


#endif // THERMPTNSAMPLER_H