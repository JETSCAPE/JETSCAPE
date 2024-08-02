
#include "ThermPtnSampler.h"
#include "JetScapeXML.h"
#include "JetScapeLogger.h"
#include "JetScapeConstants.h"
#include <vector>
#include <random>
#include <cmath>
#include <unordered_map>
#include <algorithm>
#include <gsl/gsl_sf_bessel.h>

using namespace Jetscape;

ThermalPartonSampler::ThermalPartonSampler(unsigned int ran_seed, double hydro_Tc)
    : rng_engine(ran_seed), 
      CDFTabLight(4097, std::vector<double>(2)), 
      CDFTabStrange(4097, std::vector<double>(2)),
      SetNum(false), // Set 'true' to set number of particles by hand- !!!Statistics use above temperature!!!
      SetNumLight(1000000), //If SetNum == true, this many UD quarks are generated
      SetNumStrange(0), //If SetNum == true, this many S quarks are generated
      ShuffleList(true), // Should list of particles be shuffled at the end
      degeneracy_ud(4.*6.), // Degeneracy of UD quarks
      degeneracy_s(2.*6.), // Degeneracy of S quarks
      CDFcacheAccuracy(0.003 / hbarC), // Accuracy range for cached temperature
      xmq(0.33/hbarC), // use pythia values here
      xms(0.5/hbarC), // use pythia values here
      NUMSTEP(4097), // 2^12+1, for steps of CDF Table, changes coarseness of momentum sampling
      CellDX(0.2), 
      CellDY(0.2), 
      CellDZ(0.2), 
      CellDT(0.1), 
      L(4.0), // Thickness from box edge
      W(4.0), // x/y width of box
      Time(2.0), 
      Vx(0.), // 'Uniform' no flow in x-dir
      Vy(0.), // 'Uniform' no flow in x-dir
      Vz(0.), // 'Uniform' no flow in x-dir
      hydroTc(hydro_Tc / hbarC), // converted from GeV into fm^-1
      num_ud(0), 
      num_s(0) {

    surface.clear();

    InitializeBesselK2();

	// Create the CDF table for 20 different temperatures in the range 
	// 0.95*hydroTc to 1.05*hydroTc
	for (int i = 0; i < 20; i++) {
		double T = hydroTc * (0.95 + 0.01 * i);
		createCDFCacheEntry(T, xmq, 1, CacheCDFLight, CDFTabLight);
		createCDFCacheEntry(T, xms, 2, CacheCDFStrange, CDFTabStrange);

		// Precompute part of the integral of the Fermi-Dirac distribution for 
		// light and strange quarks
		double N_eq_ud_part = FermiDiracDistributionMomentumIntegral(T, xmq);
		double N_eq_s_part = FermiDiracDistributionMomentumIntegral(T, xms);
		// Update the cache with the precomputed values
		CacheFermiDiracIntegralLight[T] = N_eq_ud_part;
		CacheFermiDiracIntegralStrange[T] = N_eq_s_part;
	}

}

void ThermalPartonSampler::createCDFCacheEntry(double& TRead, double mq, int quarkType, 
                                          std::unordered_map<double, 
										  std::vector<std::vector<double>>>& cache, 
                                          std::vector<std::vector<double>>& CDFTab) {
	// use the CDFGenerator to add the CDF data to the cache
	CDFGenerator(TRead, mq, quarkType);
	cache[TRead] = CDFTab;
}

void ThermalPartonSampler::InitializeBesselK2() {
    double K2_x_min = 0.001;
    double K2_x_max = 400.0;
    double K2_dx = 0.001;
    int K2_tb_length = static_cast<int>((K2_x_max - K2_x_min)/K2_dx) + 1;

    BesselK2.resize(K2_tb_length);
    for (int i = 0; i < K2_tb_length; i++) {
        double x = K2_x_min + i*K2_dx;
        BesselK2[i] = gsl_sf_bessel_Kn(2, x);
    }
}

double ThermalPartonSampler::BesselK2function(double arg) {
    double K2_x_min = 0.001;
    double K2_x_max = 400.0;
    double K2_dx = 0.001;
    if (arg < K2_x_min || arg > K2_x_max) {
        return gsl_sf_bessel_Kn(2, arg);
    } else {
        int idx = static_cast<int>((arg - K2_x_min)/K2_dx);
        double fraction = (arg - K2_x_min - idx*K2_dx)/K2_dx;
        return (1.0 - fraction)*BesselK2[idx] + fraction*BesselK2[idx + 1];
    }
}

double ThermalPartonSampler::FermiDiracDistributionMomentumIntegral(double T, double m) {
	int truncate_order = 10;
    double N_eq = 0.0;

    for (int n = 1; n <= truncate_order; n++) {
        double arg = n*m/T; // argument in Bessel function K_2

        // sign is (-1)^(n+1)
        double sign = (n % 2 == 0) ? -1.0 : 1.0;

        double K2 = BesselK2function(arg);
        N_eq += (sign/n)*K2;
    }
    N_eq *= m*m*T;
    return N_eq;
}

double ThermalPartonSampler::FermiDiracDistribution(double p, double m, double T){
	return 1./(exp(sqrt(p*p + m*m)/T) + 1.);
}

void ThermalPartonSampler::CDFGenerator(double T, double m, int quark) {
    double PStep;
    double PMax = 10. * T;  // CutOff for Integration
    PStep = PMax / (NUMSTEP - 1); // Stepsize in p

    std::vector<std::vector<double>>* CDFTab = nullptr;
    if (quark == 1) {
        CDFTab = &CDFTabLight;
    } else if (quark == 2) {
        CDFTab = &CDFTabStrange;
    } else {
        JSWARN << "This is not a valid quark input for CDFGenerator()";
        return;
    }

    // Initialize tabulated results for CDF
    (*CDFTab)[0][0] = 0; // For zero momentum or less...
    (*CDFTab)[0][1] = 0; // There is zero chance

    // Tabulate CDF(x) = int(0->x) PDF
    for (int i = 1; i < NUMSTEP; i++) {
        double p0 = (*CDFTab)[i - 1][0];
        double p1 = p0 + PStep;
        double Fermi0 = FermiDiracDistribution(p0, m, T);
        double Fermi1 = FermiDiracDistribution(p1, m, T);

        (*CDFTab)[i][0] = p1; // Calculate Momentum of the next step
        (*CDFTab)[i][1] = (*CDFTab)[i - 1][1] + (PStep / 2.0) * 
                                        (Fermi0 * p0 * p0 + Fermi1 * p1 * p1);
    }

    // Normalize the CDF
    double normalizationFactor = 1.0 / (*CDFTab)[NUMSTEP - 1][1];
    for (int i = 0; i < NUMSTEP; i++) {
        (*CDFTab)[i][1] *= normalizationFactor;
    }
}

std::tuple<double, double, double, double> ThermalPartonSampler::MomentumSampler(double T, int quark) {
    double PMag;
    double PMax = 10. * T;  // CutOff for Integration
    double PStep = PMax / (NUMSTEP - 1); // Stepsize in P

	// Get the appropriate CDF table from the cache closest to the target temperature
	double cachedTemp = getClosestCachedTemp((quark == 1) ? CacheCDFLight : CacheCDFStrange, T);
	// Use the temperature from the cache to get the CDF table
	std::vector<std::vector<double>>& CDFTab = (quark == 1) ? CacheCDFLight[cachedTemp] : CacheCDFStrange[cachedTemp];

	double PRoll;
	double denominator;
	while (true) {
        PRoll = ran();

        // Perform binary search to find the appropriate indices
        auto it = std::lower_bound(CDFTab.begin(), CDFTab.end(), std::vector<double>{0, PRoll},
                                   [](const std::vector<double>& a, const std::vector<double>& b) {
                                       return a[1] < b[1];
                                   });

        int floor = std::distance(CDFTab.begin(), it) - 1;
        int ceiling = floor + 1;

        if (floor < 0) floor = 0;
        if (ceiling >= NUMSTEP) ceiling = NUMSTEP - 1;

        denominator = CDFTab[ceiling][1] - CDFTab[floor][1];
        if (std::fabs(denominator) > 1e-16) {
            PMag = PStep * (PRoll - CDFTab[floor][1]) / denominator + CDFTab[floor][0];
            break;
        }
    }

    double CosT = (ran() - 0.5) * 2.0;
	double SinT = sqrt(1.0 - CosT * CosT);
    double Phi = ran() * 2.0 * pi;

    double NewX = PMag * SinT * cos(Phi);
    double NewY = PMag * SinT * sin(Phi);
    double NewZ = PMag * CosT;
    double NewP = PMag;

	return std::make_tuple(NewX, NewY, NewZ, NewP);
}

void ThermalPartonSampler::LorentzBoostMatrix(std::vector<double>& v, std::vector<std::vector<double>>& BoostMatrix, bool brick) {
	double v_square = 0.0;
	if (!brick) {
		v_square = v[1]*v[1] + v[2]*v[2] + v[3]*v[3];
		if(v_square < 10e-16){
			v_square = 10e-16;
		}
	} else {
		// This case is for the Brick, where the flow velocity can be set by hand
		v[1] = Vx;
		v[2] = Vy;
		v[3] = Vz;

		double v_square = v[1]*v[1] + v[2]*v[2] + v[3]*v[3];

		if (v_square > 1.) {
			JSWARN << "v^2 = " << v_square;
			JSWARN << "Unphysical velocity (brick flow)! Set to \"No Flow\" case";
			v[1] = 0.;
			v[2] = 0.;
			v[3] = 0.;
		}
	}
	
	v[0] = 1. / std::sqrt(1. - v_square);   // gamma - v is not four velocity

	// Lambda_u^v from (lab frame to rest frame with flow velocity)
	if (v_square < 10e-16) {
		BoostMatrix[0][0] = v[0];
		BoostMatrix[0][1] = v[0] * v[1];
		BoostMatrix[0][2] = v[0] * v[2];
		BoostMatrix[0][3] = v[0] * v[3];
		BoostMatrix[1][0] = v[0] * v[1];
		BoostMatrix[1][1] = 1.;
		BoostMatrix[1][2] = 0.;
		BoostMatrix[1][3] = 0.;
		BoostMatrix[2][0] = v[0] * v[2];
		BoostMatrix[2][1] = 0.;
		BoostMatrix[2][2] = 1.;
		BoostMatrix[2][3] = 0.;
		BoostMatrix[3][0] = v[0] * v[3];
		BoostMatrix[3][1] = 0.;
		BoostMatrix[3][2] = 0.;
		BoostMatrix[3][3] = 1.;
	} else {
		BoostMatrix[0][0] = v[0];
		BoostMatrix[0][1] = v[0] * v[1];
		BoostMatrix[0][2] = v[0] * v[2];
		BoostMatrix[0][3] = v[0] * v[3];
		BoostMatrix[1][0] = v[0] * v[1];
		BoostMatrix[1][1] = (v[0] - 1.) * v[1] * v[1] / v_square + 1.;
		BoostMatrix[1][2] = (v[0] - 1.) * v[1] * v[2] / v_square;
		BoostMatrix[1][3] = (v[0] - 1.) * v[1] * v[3] / v_square;
		BoostMatrix[2][0] = v[0] * v[2];
		BoostMatrix[2][1] = (v[0] - 1.) * v[1] * v[2] / v_square;
		BoostMatrix[2][2] = (v[0] - 1.) * v[2] * v[2] / v_square + 1.;
		BoostMatrix[2][3] = (v[0] - 1.) * v[2] * v[3] / v_square;
		BoostMatrix[3][0] = v[0] * v[3];
		BoostMatrix[3][1] = (v[0] - 1.) * v[1] * v[3] / v_square;
		BoostMatrix[3][2] = (v[0] - 1.) * v[2] * v[3] / v_square;
		BoostMatrix[3][3] = (v[0] - 1.) * v[3] * v[3] / v_square + 1.;
	}
}

std::vector<double> ThermalPartonSampler::LorentzBoost(std::vector<double>& x, std::vector<std::vector<double>>& BoostMatrix) {
	std::vector<double> result(4, 0.0);
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			result[i] += BoostMatrix[i][j] * x[j];
		}
	}
	return result;
}

void ThermalPartonSampler::SamplePartons(int Npartons, int quark, double T, bool brick, std::vector<double>& CPos, std::vector<std::vector<double>>& BoostMatrix, bool slice_boost, double eta_slice) {
	// Sample Npartons partons
	for (int i = 0; i < Npartons; i++) {
		// adding space to PList for output quarks
		std::vector<double> temp = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
		Plist.push_back(temp);

		double SpecRoll = ran(); // Probability of species die roll
		if (quark == 1) { // Sample U, Ubar, D, DBar quarks
			if (SpecRoll <= 0.25) { // UBar
				// set the last element of the vector to -2
				Plist.back()[1] = -2;
			} else if (SpecRoll <= 0.50) { // DBar
				Plist.back()[1] = -1;
			} else if (SpecRoll <= 0.75) { // D
				Plist.back()[1] = 1;
			} else { // U
				Plist.back()[1] = 2;
			}
		} else if (quark == 2) { // Sample S, Sbar quarks
			if (SpecRoll <= 0.5) { // SBar
				Plist.back()[1] = -3;
			} else { // S
				Plist.back()[1] = 3;
			}
		} else {
			JSWARN << "This is not a valid quark input for SamplePartons()";
			return;
		}

		if (!brick) {
			// Position
			// Located at x,y pos of area element
			Plist.back()[10] = CPos[0] + (ran() - 0.5) * CellDT; // Tau
			Plist.back()[7] = CPos[1] + (ran() - 0.5) * CellDX;
			Plist.back()[8] = CPos[2] + (ran() - 0.5) * CellDY;
			Plist.back()[9] = CPos[3] + (ran() - 0.5) * CellDZ;	
			if (slice_boost) { // 2+1D case
				if(std::abs(Plist.back()[9]) >= std::abs(Plist.back()[10])) {
					Plist.back()[10] = std::abs(Plist.back()[9]) + 10e-3;
				}
				double temp_t = std::sqrt(Plist.back()[10]*Plist.back()[10]-Plist.back()[9]*Plist.back()[9])
						* std::cosh(eta_slice+0.5*std::log((Plist.back()[10]+Plist.back()[9])/(Plist.back()[10]-Plist.back()[9])));
				double temp_z = std::sqrt(Plist.back()[10]*Plist.back()[10]-Plist.back()[9]*Plist.back()[9])
						* std::sinh(eta_slice+0.5*std::log((Plist.back()[10]+Plist.back()[9])/(Plist.back()[10]-Plist.back()[9])));
				Plist.back()[10] = temp_t;
				Plist.back()[9] = temp_z;
			}
		} else {
			// Position
			// Located at x,y pos of area element
			double XRoll = ran() - 0.5; // center at x=0
			double YRoll = ran() - 0.5; // center at y=0
			double ZRoll = ran() - 0.5; // center at z=0

			Plist.back()[7] = XRoll * L;
			Plist.back()[8] = YRoll * W;
			Plist.back()[9] = ZRoll * W;

			// Time
			Plist.back()[10] = Time; // Tau = L/2.: assume jet at light speed
		}

		// Sample rest frame momentum given T and mass of quark
		std::tuple<double, double, double, double> P = MomentumSampler(T, quark);
		double NewPx = std::get<0>(P);
		double NewPy = std::get<1>(P);
		double NewPz = std::get<2>(P);
		double NewP = std::get<3>(P);

		// PLab^u = g^u^t Lambda_t ^v Pres^w g_w _v  = Lambda^u _w Pres_w (with velocity -v)
		// (Lambda_u ^t with velocity v) == (Lambda^u _t with velocity -v)
		// This is boost as if particle sampled in rest frame
		double new_quark_energy;
		if (quark == 1) {
			new_quark_energy = sqrt(xmq * xmq + NewP * NewP);
		} else {
			new_quark_energy = sqrt(xms * xms + NewP * NewP);
		}
		std::vector<double> p = {new_quark_energy, NewPx, NewPy, NewPz};
		std::vector<double> momentum = LorentzBoost(p, BoostMatrix);
		Plist.back()[6] = momentum[0] * hbarC;
		Plist.back()[3] = momentum[1] * hbarC;
		Plist.back()[4] = momentum[2] * hbarC;
		Plist.back()[5] = momentum[3] * hbarC;

		// Additional information
		Plist.back()[0] = 1; // Event ID, to match jet formatting
		Plist.back()[2] = 0; // Origin, to match jet formatting
		Plist.back()[11] = 0; // Status - identifies as thermal quark
	}
}

void ThermalPartonSampler::sample_brick(){
	//preliminary parameter checks
	if(L < 0.){L = -L; JSWARN << "Negative brick length - setting to positive " << L << " fm.";}
	if(W < 0.){W = -W; JSWARN << "Negative brick width - setting to positive " << W << " fm.";}

	// Define hypersurface for brick here
	// Default t=const hypersurface
	// Vector must be covariant (negative signs on spatial components)
	std::vector<double> LFSigma(4); // LabFrame hypersurface (tau/t,x,y,eta/z=0)
	LFSigma[0] = 1.;
	LFSigma[1] = 0.;
	LFSigma[2] = 0.;
	LFSigma[3] = 0.;

	std::vector<double> CMSigma(4); // Center of mass hypersurface (tau/t,x,y,eta/z=0)
	std::vector<double> Vel(4); // Gamma & 3-velocity of cell (gamma, Vx, Vy, Vz) NOT FOUR VELOCITY

	// Calculated quantities in cells
	std::vector<std::vector<double>> LorBoost(4, std::vector<double>(4)); // Lorentz boost defined as used - form is always Lambda_u^v
	LorentzBoostMatrix(Vel, LorBoost, true);

	// Lambda_u^v Sigma_v = CMSigma_u
	// Caution: CMSigma is a covariant vector
	std::vector<double> SigmaBoosted = LorentzBoost(LFSigma, LorBoost);
	CMSigma[0] = SigmaBoosted[0];
	CMSigma[1] = SigmaBoosted[1];
	CMSigma[2] = SigmaBoosted[2];
	CMSigma[3] = SigmaBoosted[3];

	// GAUSSIAN INTEGRALS <n> = int f(p)d3p
	double dSigma_dot_u = CMSigma[0] * Vel[0] - CMSigma[1] * Vel[1] - CMSigma[2] * Vel[2] - CMSigma[3] * Vel[3];

	// get the closest cached temperature to the hydrodynamic temperature
	double TCachedLight = getClosestCachedTemp(CacheFermiDiracIntegralLight, hydroTc);
	double TCachedStrange = getClosestCachedTemp(CacheFermiDiracIntegralStrange, hydroTc);

	// get the precomputed part of the Fermi-Dirac integrals for light and strange quarks from cache
	// U, D, UBAR, DBAR QUARKS
	// <N> = V <n>
	double NumUD = (CacheFermiDiracIntegralLight[TCachedLight] * degeneracy_ud * dSigma_dot_u 
				/ (2.*pi*pi)) * L * W * W;

	// S, SBAR QUARKS
	// <N> = V <n>
	double NumS = (CacheFermiDiracIntegralStrange[TCachedStrange] * degeneracy_s * dSigma_dot_u 
				/ (2.*pi*pi)) * L * W * W;

	std::poisson_distribution<int> poisson_ud(NumUD);
	std::poisson_distribution<int> poisson_s(NumS);

	// In case of overwriting
	if(SetNum){
		NumUD = SetNumLight;
		NumS = SetNumStrange;
	}

	// Generating light quarks
	int GeneratedParticles = poisson_ud(getRandomGenerator());
	num_ud += GeneratedParticles;
	// dummy array of length 4 for CPos
	std::vector<double> CPos(4, 0.0);
	SamplePartons(GeneratedParticles, 1, hydroTc, true, CPos, LorBoost, false, 0.);

	// Generate strange quarks
	GeneratedParticles = poisson_s(getRandomGenerator());
	num_s += GeneratedParticles;
	SamplePartons(GeneratedParticles, 2, hydroTc, true, CPos, LorBoost, false, 0.);

	JSDEBUG << "Light particles: " << num_ud;
	JSDEBUG << "Strange particles: " << num_s;

	// Shuffling PList
	if (!Plist.empty() && ShuffleList) {
    	std::shuffle(Plist.begin(), Plist.end(), getRandomGenerator());
	}

	//print Plist for testing
	/*std::cout << std::setprecision(5);
  	std::ofstream thermalP;
  	thermalP.open("thermal_partons.dat", std::ios::app);
	for(int p=0; p < Plist.size(); p++){
		thermalP << Plist[p][0] << " " << Plist[p][1] << " " << Plist[p][2]
		<< " " << Plist[p][3] << " " << Plist[p][4] << " " << Plist[p][5]
		<< " " << Plist[p][6] << " " << Plist[p][7] << " " << Plist[p][8]
		<< " " << Plist[p][9] << " " << Plist[p][10] << " " << Plist[p][11]
		<< "\n";
	}
	thermalP.close();*/
}

void ThermalPartonSampler::sample_3p1d(bool Cartesian_hydro){

	// Loop over the freeze-out surface
	for(int iS=0; iS<surface.size(); ++iS){
		std::vector<double> CPos(4); // Position of the current cell (tau/t, x, y , eta/z=0)
		double tau_pos = surface[iS][0]; // proper time from position of cell
		CPos[1] = surface[iS][1];
		CPos[2] = surface[iS][2];
		double eta_pos = surface[iS][3]; // eta from position of cell, we need t,x,y,z
		
		//this is also tau, x, y, eta
		std::vector<double> LFSigma(4); // LabFrame hypersurface (tau/t,x,y,eta/z), expect Sigma_mu
		std::vector<double> CMSigma(4); // Center of mass hypersurface (tau/t,x,y,eta/z)
		double tau_sur = surface[iS][4]; // proper time from normal vector of surface
		LFSigma[1] = surface[iS][5];
		LFSigma[2] = surface[iS][6];
		double eta_sur = surface[iS][7]; // eta from normal vector of surface
		
		double TRead = surface[iS][8] / hbarC; // Temperature of cell -> convert to fm^-1
		
		std::vector<double> Vel(4); // Gamma & 3-Velocity of cell (gamma, Vx, Vy, Vz) NOT FOUR VELOCITY
		Vel[1] = surface[iS][9];
		Vel[2] = surface[iS][10];
		Vel[3] = surface[iS][11];

		if (Cartesian_hydro == false){
			//getting t,z from tau, eta
			CPos[0] = tau_pos*std::cosh(eta_pos);
			CPos[3] = tau_pos*std::sinh(eta_pos);
			//transform surface vector to Minkowski coordinates
			double cosh_eta_pos = std::cosh(eta_pos);
			double sinh_eta_pos = std::sinh(eta_pos);
			LFSigma[0] = cosh_eta_pos*tau_sur - (sinh_eta_pos / tau_pos) * eta_sur;
			LFSigma[3] = -sinh_eta_pos*tau_sur + (cosh_eta_pos / tau_pos) * eta_sur;
		}else{
			CPos[0] = tau_pos;
			CPos[3] = eta_pos;
			LFSigma[0] = tau_sur;
			LFSigma[3] = eta_sur;
		}

		std::vector<std::vector<double>> LorBoost(4, std::vector<double>(4)); // Lorentz boost defined as used - Form is always Lambda_u^v
		LorentzBoostMatrix(Vel, LorBoost, false);

		// Lambda_u^v Sigma_v = CMSigma_u
		std::vector<double> SigmaBoosted = LorentzBoost(LFSigma, LorBoost);
		CMSigma[0] = SigmaBoosted[0];
		CMSigma[1] = SigmaBoosted[1];
		CMSigma[2] = SigmaBoosted[2];
		CMSigma[3] = SigmaBoosted[3];

		// INTEGRAL <n> = int f(p)d3p
		// get the precomputed part of the Fermi-Dirac integrals for light and strange quarks from cache
		double dSigma_dot_u = CMSigma[0] * Vel[0] - CMSigma[1] * Vel[1] - CMSigma[2] * Vel[2] - CMSigma[3] * Vel[3];
		// get the precomputed part of the Fermi-Dirac integrals for light and strange quarks from cache closest to the temperature of the cell
		double TCacheLight = getClosestCachedTemp(CacheFermiDiracIntegralLight, TRead);
		double TCacheStrange = getClosestCachedTemp(CacheFermiDiracIntegralStrange, TRead);
		double NumLight = CacheFermiDiracIntegralLight[TCacheLight] * degeneracy_ud * dSigma_dot_u / (2.*pi*pi);
		double NumStrange = CacheFermiDiracIntegralStrange[TCacheStrange] * degeneracy_s * dSigma_dot_u / (2.*pi*pi);

		// Generating light quarks
		std::poisson_distribution<int> poisson_ud(NumLight);
		int GeneratedParticles_ud = poisson_ud(getRandomGenerator()); // Initialize particles created in this cell
		SamplePartons(GeneratedParticles_ud, 1, TRead, false, CPos, LorBoost, false, 0.);
		num_ud += GeneratedParticles_ud;

		// Generate s quarks
		std::poisson_distribution<int> poisson_s(NumStrange);
		int GeneratedParticles_s = poisson_s(getRandomGenerator()); //Initialize particles created in this cell
		SamplePartons(GeneratedParticles_s, 2, TRead, false, CPos, LorBoost, false, 0.);
		num_s += GeneratedParticles_s;
	}

	JSDEBUG << "Light particles: " << num_ud;
	JSDEBUG << "Strange particles: " << num_s;

	//Shuffling PList
	if (!Plist.empty() && ShuffleList) {
    	std::shuffle(Plist.begin(), Plist.end(), getRandomGenerator());
	}

	//print Plist for testing
	/*std::cout << std::setprecision(5);
  	std::ofstream thermalP;
  	thermalP.open("thermal_partons.dat", std::ios::app);
	for(int p=0; p < Plist.size(); p++){
		thermalP << Plist[p][0] << " " << Plist[p][1] << " " << Plist[p][2]
		<< " " << Plist[p][3] << " " << Plist[p][4] << " " << Plist[p][5]
		<< " " << Plist[p][6] << " " << Plist[p][7] << " " << Plist[p][8]
		<< " " << Plist[p][9] << " " << Plist[p][10] << " " << Plist[p][11]
		<< "\n";
	}
	thermalP.close();*/
}


void ThermalPartonSampler::sample_2p1d(double eta_max){

	// precompute CDFs and store them in a cache
	std::unordered_map<double, std::vector<std::vector<double>>> cdfLightCache;
	std::unordered_map<double, std::vector<std::vector<double>>> cdfStrangeCache;

	std::vector<double> NumberLightQuarks;
	std::vector<double> NumberStrangeQuarks;

	double d_eta = CellDZ;
	int N_slices = std::floor((eta_max / d_eta)-0.5)+1;
	for(int slice=1; slice <= (2*N_slices+1); slice++){
		double eta_slice = (slice-N_slices-1)*d_eta;

		JSINFO << "Sampling thermal partons for pseudorapidity slice " << slice
		<< " (eta_min = " << eta_slice-(d_eta/2.) << ", eta_max = "
		<< eta_slice+(d_eta/2.) << ")";

		for(int iS=0; iS<surface.size(); ++iS){
			std::vector<double> CPos(4); // Position of the current cell (tau/t, x, y , eta/z=0)
			double tau_pos = surface[iS][0]; // proper time from position of cell
			CPos[1] = surface[iS][1];
			CPos[2] = surface[iS][2];
			double eta_pos = surface[iS][3]; // eta from position of cell, we need t,x,y,z
			
			//this is also tau, x, y, eta
			std::vector<double> LFSigma(4); // LabFrame hypersurface (tau/t,x,y,eta/z=0), expect Sigma_mu
			std::vector<double> CMSigma(4); // Center of mass hypersurface (tau/t,x,y,eta/z=0)
			double tau_sur = surface[iS][4]; // proper time from normal vector of surface
			LFSigma[1] = surface[iS][5];
			LFSigma[2] = surface[iS][6];
			double eta_sur = surface[iS][7]; // eta from normal vector of surface
			
			double TRead = surface[iS][8] / hbarC; // Temperature of cell -> convert to fm^-1
			
			std::vector<double> Vel(4); // Gamma & 3-Velocity of cell (gamma, Vx, Vy, Vz) NOT FOUR VELOCITY
			Vel[1] = surface[iS][9];
			Vel[2] = surface[iS][10];
			Vel[3] = surface[iS][11];

			//getting t,z from tau, eta
			CPos[0] = tau_pos*std::cosh(eta_pos);
			CPos[3] = tau_pos*std::sinh(eta_pos);
			//transform surface vector to Minkowski coordinates
			double cosh_eta_pos = std::cosh(eta_pos);
			double sinh_eta_pos = std::sinh(eta_pos);
			LFSigma[0] = cosh_eta_pos*tau_sur - (sinh_eta_pos / tau_pos) * eta_sur;
			LFSigma[3] = -sinh_eta_pos*tau_sur + (cosh_eta_pos / tau_pos) * eta_sur;

			CellDZ = CPos[0] * 2. * std::sinh(d_eta/2.);

			std::vector<std::vector<double>> LorBoost(4, std::vector<double>(4)); // Lorentz boost defined as used - Form is always Lambda_u^v
			LorentzBoostMatrix(Vel, LorBoost, false);

			// Lambda_u^v Sigma_v = CMSigma_u
			std::vector<double> SigmaBoosted = LorentzBoost(LFSigma, LorBoost);
			CMSigma[0] = SigmaBoosted[0];
			CMSigma[1] = SigmaBoosted[1];
			CMSigma[2] = SigmaBoosted[2];
			CMSigma[3] = SigmaBoosted[3];

			// INTEGRAL <n> = int f(p)d3p
			// get the precomputed part of the Fermi-Dirac integrals for light and strange quarks from cache
			double dSigma_dot_u = CMSigma[0] * Vel[0] - CMSigma[1] * Vel[1] - CMSigma[2] * Vel[2] - CMSigma[3] * Vel[3];
			// get the precomputed part of the Fermi-Dirac integrals for light and strange quarks from cache closest to the temperature of the cell
			double TCacheLight = getClosestCachedTemp(CacheFermiDiracIntegralLight, TRead);
			double TCacheStrange = getClosestCachedTemp(CacheFermiDiracIntegralStrange, TRead);

			double NumLight = CacheFermiDiracIntegralLight[TCacheLight] * degeneracy_ud * dSigma_dot_u / (2.*pi*pi);
			double NumStrange = CacheFermiDiracIntegralStrange[TCacheStrange] * degeneracy_s * dSigma_dot_u / (2.*pi*pi);

			// Velocity for the boost to the slice
			Vel[1] /= cosh(eta_slice);
			Vel[2] /= cosh(eta_slice);
			Vel[3] = tanh(eta_slice);
			// Update the Lorentz boost matrix for the slice
			LorentzBoostMatrix(Vel, LorBoost, false);

			// Generating light quarks
			std::poisson_distribution<int> poisson_ud(NumLight);
			int GeneratedParticles_ud = poisson_ud(getRandomGenerator()); // Initialize particles created in this cell
			SamplePartons(GeneratedParticles_ud, 1, TRead, false, CPos, LorBoost, true, eta_slice);
			num_ud += GeneratedParticles_ud;

			// Generate s quarks
			std::poisson_distribution<int> poisson_s(NumStrange);
			int GeneratedParticles_s = poisson_s(getRandomGenerator()); //Initialize particles created in this cell
			SamplePartons(GeneratedParticles_s, 2, TRead, false, CPos, LorBoost, true, eta_slice);
			num_s += GeneratedParticles_s;
		}

	}

	JSDEBUG << "Light particles: " << num_ud;
	JSDEBUG << "Strange particles: " << num_s;

	if (!Plist.empty() && ShuffleList) {
    	std::shuffle(Plist.begin(), Plist.end(), getRandomGenerator());
	}

	//print Plist for testing
	/*std::cout << std::setprecision(5);
  	std::ofstream thermalP;
  	thermalP.open("thermal_partons.dat", std::ios::app);
	for(int p=0; p < Plist.size(); p++){
		thermalP << Plist[p][0] << " " << Plist[p][1] << " " << Plist[p][2]
		<< " " << Plist[p][3] << " " << Plist[p][4] << " " << Plist[p][5]
		<< " " << Plist[p][6] << " " << Plist[p][7] << " " << Plist[p][8]
		<< " " << Plist[p][9] << " " << Plist[p][10] << " " << Plist[p][11]
		<< "\n";
	}
	thermalP.close();*/
}
