#ifndef THERMPTNSAMPLER_H
#define THERMPTNSAMPLER_H

#include "JetScapeLogger.h"
#include <vector>
#include <random>


using namespace Jetscape;

class ThermalPartonSampler
{
 public:
	/**
     * @brief Constructor that initializes the random number generator 
	 * and thermal freeze-out temperature.
     * 
     * @param ran_seed Random seed for the random number generator.
     * @param hydro_Tc Critical temperature for hydrodynamic evolution (in GeV).
     */
	ThermalPartonSampler(unsigned int ran_seed, double hydro_Tc);

	/**
	 * @brief Samples partons within a brick-like hydrodynamics setup.
	 * 
	 * Generates partons according to the specified parameters within a brick-shaped 
	 * hypersurface defined by the length (L) and width (W) parameters. Partons are 
	 * sampled according to Fermi-Dirac distributions for light (u, d) and strange 
	 * quarks.
	 * 
	 * @details
	 * - Checks and corrects negative length and width parameters.
	 * - Defines the hypersurface vectors in the lab frame and boosts them to the 
	 *   rest frame.
	 * - Calculates the number of quarks (light and strange) using Fermi-Dirac 
	 *   distributions and volume.
	 * - Generates partons for light (u, d) and strange quarks based on Poisson 
	 *   distributions of expected numbers.
	 * - Adjusts the number of generated particles if overridden by `SetNumLight` 
	 *   and `SetNumStrange`.
	 * - Shuffles the generated `Plist` if `ShuffleList` is true.
	 */
	void sample_brick();
	
	/**
	 * @brief Samples partons on a freeze-out hypersurface in 2+1D setup.
	 * 
	 * Samples light (u, d, u-bar, d-bar) and strange (s, s-bar) quarks based on 
	 * local temperature and velocity. CDFs are precomputed and cached to optimize 
	 * sampling.
	 * 
	 * @param eta_max Maximum pseudorapidity to sample over.
	 */
	void sample_2p1d(double eta_max);

	/**
	 * @brief Samples partons on a freeze-out hypersurface in 3+1D setup.
	 * 
	 * Samples light (u, d, u-bar, d-bar) and strange (s, s-bar) quarks based on 
	 * local temperature and velocity. CDFs are precomputed and cached to optimize 
	 * sampling.
	 * 
	 * @param Cartesian_hydro Boolean indicating whether Cartesian hydrodynamic 
	 * coordinates are used. If false, uses Milne coordinates.
	 */
	void sample_3p1d(bool Cartesian_hydro=false);

	/**
     * @brief Sets the hypersurface data.
     * 
     * @param surf_in 2D vector representing the hypersurface data.
     */
	void set_hypersurface(std::vector<std::vector<double>> surf_in){surface = surf_in;}

	/**
     * @brief Sets the brick dimensions.
     * 
     * @param len_bri Length of the brick in fm.
     * @param wid_bri Width of the brick in fm.
     */
	void brick_length_width(double len_bri, double wid_bri){L = 2.*len_bri + 4.; W = 2.*wid_bri + 4.; Time = len_bri;} // +4 gives 2fm additional brick in each direction
	
	/**
     * @brief Sets the flow velocity of the brick.
     * 
     * @param vx_in Flow velocity component in the x-direction.
     * @param vy_in Flow velocity component in the y-direction.
     * @param vz_in Flow velocity component in the z-direction.
     */
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

	// getter functions for the brick parameters
	double getL(){return L;}
	double getW(){return W;}
	double getTime(){return Time;}
	double getVx(){return Vx;}
	double getVy(){return Vy;}
	double getVz(){return Vz;}

	// getter functions for the cached CDFs for testing
	std::unordered_map<double, std::vector<std::vector<double>>>& getCacheCDFLight(){return CacheCDFLight;}
	std::unordered_map<double, std::vector<std::vector<double>>>& getCacheCDFStrange(){return CacheCDFStrange;}

	/**
	 * @brief Computes the modified Bessel function of the second kind, K2(x), 
	 * using a lookup table for efficiency.
	 * 
	 * This function returns the value of the modified Bessel function K2(x) for a 
	 * given argument. If the argument is within the precomputed range, it uses 
	 * linear interpolation between the nearest values in the lookup table to 
	 * provide a fast approximation. If the argument is outside the precomputed 
	 * range, it directly computes the value using the GSL library.
	 * 
	 * @param arg The argument for which the Bessel function value is to be 
	 * computed.
	 * @return The value of K2(arg).
	 */
	double BesselK2function(double arg);

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
	
	/**
     * @brief Generates a random number between 0 and 1.
     * 
     * @return A random double in the range [0.0, 1.0].
     */
    double ran() {return distribution(rng_engine);}
	
	/**
     * @brief Gets the random number generator.
     * 
     * @return A reference to the Mersenne Twister random number generator.
     */
    std::mt19937_64& getRandomGenerator() {return rng_engine;}

	
	// HyperSurface
	std::vector<std::vector<double>> surface;

	// Vector to pre-tabulate the Bessel function K2
	std::vector<double> BesselK2;

	/**
	 * @brief Initializes the lookup table for the modified Bessel function of the 
	 * second kind, K2(x).
	 * 
	 * This function calculates K2(x) for x values from 0.001 to 400.0 with an 
	 * increment of 0.001, and stores the results in the BesselK2 vector for 
	 * efficient lookup.
	 */
	void InitializeBesselK2();

	/**
	 * @brief Computes the integral of the Fermi-Dirac distribution over momentum.
	 * 
	 * This function calculates the integral of the Fermi-Dirac distribution, which 
	 * is used to determine the equilibrium number density of fermions at a given 
	 * temperature (T) and mass (m). The calculation is performed using a truncated 
	 * series expansion.
	 * 
	 * @param T The temperature.
	 * @param m The mass of the fermion.
	 * @return The computed equilibrium number density, N_eq.
	 */
	double FermiDiracDistributionMomentumIntegral(double T, double m);
	
	/**
	 * @brief Computes the Fermi-Dirac distribution function.
	 * 
	 * This function calculates the value of the Fermi-Dirac distribution for a 
	 * given momentum (p), mass (m), and temperature (T).
	 * The Fermi-Dirac distribution describes the probability of occupancy of a 
	 * quantum state by a fermion at thermal equilibrium.
	 * 
	 * @param p The momentum of the fermion.
	 * @param m The mass of the fermion.
	 * @param T The temperature.
	 * @return The value of the Fermi-Dirac distribution function.
	 */
	double FermiDiracDistribution(double p, double m, double T);

	/**
	 * @brief Generates a cumulative distribution function (CDF) table for a given 
	 * temperature, mass, and quark type.
	 * 
	 * This function computes the cumulative distribution function (CDF) for the 
	 * Fermi-Dirac distribution at a specified temperature (T) and mass (m) for 
	 * either light or strange quarks. The CDF values are stored in the appropriate 
	 * table (CDFTabLight or CDFTabStrange).
	 * 
	 * @param T The temperature.
	 * @param m The mass of the quark.
	 * @param quark An integer indicating the type of quark (1 for light quark, 
	 * 2 for strange quark).
	 */
	void CDFGenerator(double T, double m, int quark);
	
	/**
	 * @brief Samples a momentum vector from the thermal distribution for a 
	 * specified quark type.
	 * 
	 * This function generates a random momentum vector for a quark (light or 
	 * strange) using the cumulative distribution function (CDF) table. The sampling 
	 * is performed by generating a random probability and finding the corresponding 
	 * momentum magnitude using binary search. The momentum vector is then generated
	 * with a random direction.
	 * 
	 * @param T The temperature in fm^{-1}.
	 * @param quark An integer indicating the type of quark (1 for light quark, 
	 * otherwise for strange quark).
	 * @return A tuple containing the components of the momentum vector 
	 * (NewX, NewY, NewZ) and the magnitude (NewP).
	 */
	std::tuple<double, double, double, double> MomentumSampler(double T, int quark);

	/**
	 * @brief Computes the Lorentz boost matrix corresponding to a given velocity 
	 * vector.
	 * 
	 * This function calculates the Lorentz boost matrix for transforming from the 
	 * laboratory frame to the rest frame with a specified flow velocity vector. 
	 * 
	 * @param v A vector containing the components of the velocity (v0, vx, vy, vz).
	 * @param BoostMatrix A 4x4 matrix that will be filled with the Lorentz boost 
	 * matrix.
	 * @param brick A boolean flag indicating whether the velocity vector is for 
	 * a "brick" flow (true) or a general case (false).
	 */
    void LorentzBoostMatrix(std::vector<double>& v, std::vector<std::vector<double>>& BoostMatrix, bool brick);

    /**
	 * @brief Applies a Lorentz boost transformation to a 4-vector using a given
	 * boost matrix.
	 * 
	 * This function performs a Lorentz boost transformation on a 4-vector `x` 
	 * using the specified 4x4 boost matrix `BoostMatrix`. 
	 * 
	 * @param x A 4-element vector representing the original 4-vector to be boosted.
	 * @param BoostMatrix A 4x4 matrix representing the Lorentz boost matrix.
	 * @return A 4-element vector containing the components of the boosted 4-vector.
	 */
    std::vector<double> LorentzBoost(std::vector<double>& x, std::vector<std::vector<double>>& BoostMatrix);

    /**
	 * @brief Samples partons according to specified parameters.
	 * 
	 * Generates a specified number of partons (`Npartons`) with properties 
	 * determined by the temperature (`T`), quark type (`quark`), and configuration 
	 * (`brick`). Positions and momenta are stored in the member variable `Plist`, 
	 * which is a vector of vectors representing each parton's properties.
	 * 
	 * @param Npartons Number of partons to sample.
	 * @param quark Type of quark to sample (1 for U/D quarks, 2 for S quarks).
	 * @param T Temperature of the system.
	 * @param brick Flag indicating whether to use a brick configuration.
	 * @param CPos Center position in spatial coordinates.
	 * @param BoostMatrix Lorentz boost matrix for momentum transformation.
	 * @param slice_boost Flag indicating 2+1D boost configuration.
	 * @param eta_slice Pseudo-rapidity for 2+1D boost.
	 */
    void SamplePartons(int Npartons, int quark, double T, bool brick, std::vector<double>& CPos, std::vector<std::vector<double>>& BoostMatrix, bool slice_boost, double eta_slice);

	/**
	 * @brief Creates an entry of the cumulative distribution function (CDF) data 
	 * in a cache.
	 * 
	 * @param TRead Input temperature to read or update in the cache.
	 * @param mq Quark mass associated with the CDF data.
	 * @param quarkType Type of quark (1 for light, 2 for strange).
	 * @param cache Reference to the cache storing CDF data for different 
	 * temperatures.
	 * @param CDFTab Reference to the vector storing the CDF table for the current 
	 * temperature.
	 */
	void createCDFCacheEntry(double& TRead, double xmq, int quarkType, 
                        std::unordered_map<double,
						std::vector<std::vector<double>>>& cache, 
                        std::vector<std::vector<double>>& CDFTab);

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

/**
 * @brief Finds the temperature closest to a target temperature in a cached map
 * for the CDF tables.
 * 
 * Searches through the provided cache of temperatures (`cache`) to determine
 * the temperature that is closest to the specified `targetTemp`.
 * 
 * @param cache A map where keys are cached temperatures and values are 
 * associated data.
 * @param targetTemp The target temperature for which the closest cached 
 * temperature is sought.
 * @return The cached temperature that is closest to `targetTemp`, or -1.0 if 
 * the cache is empty.
 */
template<typename Cache>
double getClosestCachedTemp(const Cache& cache, double targetTemp) {
    double closestTemp = -1.0;
    double minDistance = std::numeric_limits<double>::max();
    for (const auto& entry : cache) {
        double cachedTemp = entry.first;
        double distance = std::fabs(targetTemp - cachedTemp);
        if (distance < minDistance) {
            minDistance = distance;
            closestTemp = cachedTemp;
        }
    }
    return closestTemp;
}


#endif // THERMPTNSAMPLER_H