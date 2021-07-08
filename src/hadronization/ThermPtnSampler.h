#ifndef THERMPTNSAMPLER_H
#define THERMPTNSAMPLER_H

#include "JetScapeLogger.h"
#include <vector>
#include <random>


using namespace Jetscape;

class thermptnsampler
{
 public:
	//constructor - initializes variables
	thermptnsampler();
	
	//sampler - samples thermal partons (u,d,s) from brick or 2+1d / 3+1d hydro
	void samplebrick(); void sample_2p1d(); void sample_3p1d();
	
	//setters for params
	void brick_LW( double len_bri, double wid_bri){                  L = len_bri; W = wid_bri; offset = L;}
	void brick_LWo(double len_bri, double wid_bri, double offset_in){L = len_bri; W = wid_bri; offset = offset_in;}
	void brick_off(double offset_in){offset = offset_in;}
	void brick_flow(double vx_in, double vy_in, double vz_in){Vx = vx_in; Vy = vy_in; Vz = vz_in;}
	void brick_seed(unsigned int seed_in){ran_seed = seed_in;}
	
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
	int th_nL(){return num_ud;} int th_nS(){return num_s;}
	
	
 private:
	
	//thermal parton list
	//double Plist[1000000][12];
	std::vector<std::vector<double>> Plist;	// List of Particles ( [0]->Event Number; [1]->Particle ID; [2]->Origin; [3-6]->Px,Py,Pz,E; [7-10]->x,y,z,t; [11]->Particle Status )
											// Same format as in shower data, Event number is always 1, Origin is always 0, particle statis is always 0 (Indicates a thermal quark)
	
	/* Functions*/
	bool SplitSort (double goal, int floor, int ceiling, int quark);	
	void MCSampler(double Temp, double M, int quark);					//Samples momentum distribution, redefines NewP, NewX, NewY, NewZ variables- momentum in rest frame
	double FermiPDF (double Pc, double Mc, double Tc, double muc);		//Gives form of Fermi-Dirac distribution (no normalization factors!). Mu is currently given by 0 everywhere
	//double BosePDF (double P, double M, double T, double mu);			
	//double BoltzPDF (double P, double M, double T, double mu);		
	void CDFGen(double Temp, double M, int quark);						//generates Cumulative distribution in thermal rest frame "quark" is 1 for light and 2 for strange
	double ran(); unsigned int ran_seed; std::mt19937_64 eng;			//RNG - Mersenne Twist - 64 bit
	
	/* Static Parameters, Do not change */
	double PIE, GEVFM, muPi0;
	
	/* Adjustable Parameters */
	double ML, MS, /*DEta,*/ T, NUMSTEP;
	double CellDX, CellDY, CellDZ, CellDT;
	
	/* Gaussian Weights and Abscissa (50 pt) */
	/* Used in numeric integration */
	int GPoints = 50;		// Amount of Gaussian Points for Quadrature
	double GWeight[50];		/* Gaussian Weights for Integration */
	double GAbs[50];		/* Gaussian Abscissas for Integration */
	
	/* Flags */
	bool SetNum, SetNumLight, SetNumStrange, ShuffleList;
	
	/* Global Vars between functions */
	double NewX, NewY, NewZ, NewP;
	int num_ud, num_s;
	//double CDFTabLight[NUMSTEP][2];   // Cumulative Distribution for thermal Light Quarks
	//double CDFTabStrange[NUMSTEP][2]; // Cumulative Distribution for Squarks
	std::vector<std::vector<double>> CDFTabLight;	// Cumulative Distribution for thermal Light Quarks
	std::vector<std::vector<double>> CDFTabStrange;	// Cumulative Distribution for Squarks
	
	/* Brick Info */
	//L is the length of the brick sampled for thermal partons
	//W is the width of the face of the brick - should be scaled somewhat against L
	//offset is where the leading edge of the brick is (eg. a L=4fm brick can be from 0-4fm(offset=4fm) or from 2-6fm(offset=6fm))
	//Vx,Vy,Vz is a flow given to the brick
	double L, W, offset, Vx, Vy, Vz;
	
};


#endif // THERMPTNSAMPLER_H